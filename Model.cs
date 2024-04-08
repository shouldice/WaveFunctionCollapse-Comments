// Copyright (C) 2016 Maxim Gumin, The MIT License (MIT)

using System;

abstract class Model
{
    protected bool[][] outputArray;

    // Dimensions: [Direction], [Pattern], [Pattern]
    // For each direction, for each pattern, a list of "agreeing" patterns
    // AKA "propagator"
    protected int[][][] matchingPatternsInDir;
    
    // Dimensions: [Output Position], [Patterns], [Direction].
    // For each output position, for each pattern, an array of 4 directions.
    // It's filled with this:
    //          compatible[i][t][d] = propagator[opposite[d]][t].Length;
    // ... where opposite is used to 180º a direction.
    // So compatible's direction ints are initialized to the NUMBER of agreeing patterns in the oppoosite direction?
    // AKA "compatible"
    int[][][] matchingPatternsCount;
    
    protected int[] observed;

    (int, int)[] stack;
    int stackPointer, observedSoFar;

    protected int outputWidth, outputHeight, patternCount, PatternSize;
    protected bool periodic, ground;

    protected double[] weights;
    double[] weightLogWeights, distribution;

    protected int[] sumsOfOnes;
    double sumOfWeights, sumOfWeightLogWeights, startingEntropy;
    protected double[] sumsOfWeights, sumsOfWeightLogWeights, entropies;

    public enum Heuristic { Entropy, MRV, Scanline };
    Heuristic heuristic;

    protected Model(int width, int height, int patternSize, bool periodic, Heuristic heuristic)
    {
        outputWidth = width;
        outputHeight = height;
        this.PatternSize = patternSize;
        this.periodic = periodic;
        this.heuristic = heuristic;
    }

    void Init()
    {
        outputArray = new bool[outputWidth * outputHeight][];
        matchingPatternsCount = new int[outputArray.Length][][];
        for (int i = 0; i < outputArray.Length; i++)
        {
            outputArray[i] = new bool[patternCount];
            matchingPatternsCount[i] = new int[patternCount][];
            for (int t = 0; t < patternCount; t++) matchingPatternsCount[i][t] = new int[4];
        }
        distribution = new double[patternCount];
        observed = new int[outputWidth * outputHeight];

        weightLogWeights = new double[patternCount];
        sumOfWeights = 0;
        sumOfWeightLogWeights = 0;

        for (int t = 0; t < patternCount; t++)
        {
            weightLogWeights[t] = weights[t] * Math.Log(weights[t]);
            sumOfWeights += weights[t];
            sumOfWeightLogWeights += weightLogWeights[t];
        }

        startingEntropy = Math.Log(sumOfWeights) - sumOfWeightLogWeights / sumOfWeights;

        sumsOfOnes = new int[outputWidth * outputHeight];
        sumsOfWeights = new double[outputWidth * outputHeight];
        sumsOfWeightLogWeights = new double[outputWidth * outputHeight];
        entropies = new double[outputWidth * outputHeight];

        stack = new (int, int)[outputArray.Length * patternCount];
        stackPointer = 0;
    }

    public bool Run(int seed, int limit)
    {
        if (outputArray == null) Init();

        Clear();
        Random random = new(seed);

        for (int l = 0; l < limit || limit < 0; l++)
        {
            int node = NextUnobservedNode(random);
            if (node >= 0)
            {
                Observe(node, random);
                bool success = Propagate();
                if (!success) return false;
            }
            else
            {
                for (int i = 0; i < outputArray.Length; i++) for (int t = 0; t < patternCount; t++) if (outputArray[i][t]) { observed[i] = t; break; }
                return true;
            }
        }

        return true;
    }

    int NextUnobservedNode(Random random)
    {
        if (heuristic == Heuristic.Scanline)
        {
            for (int i = observedSoFar; i < outputArray.Length; i++)
            {
                if (!periodic && (i % outputWidth + PatternSize > outputWidth || i / outputWidth + PatternSize > outputHeight)) continue;
                if (sumsOfOnes[i] > 1)
                {
                    observedSoFar = i + 1;
                    return i;
                }
            }
            return -1;
        }

        double min = 1E+4;
        int argmin = -1;
        for (int i = 0; i < outputArray.Length; i++)
        {
            if (!periodic && (i % outputWidth + PatternSize > outputWidth || i / outputWidth + PatternSize > outputHeight)) continue;
            int remainingValues = sumsOfOnes[i];
            double entropy = heuristic == Heuristic.Entropy ? entropies[i] : remainingValues;
            if (remainingValues > 1 && entropy <= min)
            {
                double noise = 1E-6 * random.NextDouble();
                if (entropy + noise < min)
                {
                    min = entropy + noise;
                    argmin = i;
                }
            }
        }
        return argmin;
    }

    void Observe(int node, Random random)
    {
        bool[] w = outputArray[node];
        for (int t = 0; t < patternCount; t++) distribution[t] = w[t] ? weights[t] : 0.0;
        int r = distribution.Random(random.NextDouble());
        for (int t = 0; t < patternCount; t++) if (w[t] != (t == r)) Ban(node, t);
    }

    bool Propagate()
    {
        while (stackPointer > 0)
        {
            (int location, int pattern) = stack[stackPointer - 1];
            stackPointer--;

            int x1 = location % outputWidth;
            int y1 = location / outputWidth;

            // d is an index into orthogonal directions?
            for (int direction = 0; direction < 4; direction++)
            {
                // x1,y1 is the source? x2,y2 is the adjacent spot?
                int x2 = x1 + dx[direction];
                int y2 = y1 + dy[direction];
                if (!periodic && (x2 < 0 || y2 < 0 || x2 + PatternSize > outputWidth || y2 + PatternSize > outputHeight)) continue;

                // wrapping
                if (x2 < 0) x2 += outputWidth;
                else if (x2 >= outputWidth) x2 -= outputWidth;
                if (y2 < 0) y2 += outputHeight;
                else if (y2 >= outputHeight) y2 -= outputHeight;

                // generate an unfurled index (i2)
                int i2 = x2 + y2 * outputWidth;
                
                // matchingPatternsInDir/"propogator" is a 3D array: for each direction, for each pattern, a list of matching patterns
                // matchingPatternsList = a list of matching patterns...........
                int[] matchingPatternList = matchingPatternsInDir[direction][pattern];
                
                // "compatible"/matchingPatternsCount is another 2D array: get the 2D array at i2 (the location of the adjacent column)
                // p = [patterns],[direction] a count of how many matching patterns in that direction.
                int[][] p = matchingPatternsCount[i2];

                for (int l = 0; l < matchingPatternList.Length; l++)
                {
                    int t2 = matchingPatternList[l];
                    int[] comp = p[t2];

                    comp[direction]--;
                    if (comp[direction] == 0) Ban(i2, t2);
                }
            }
        }

        return sumsOfOnes[0] > 0;
    }

    // "Ban" as in disallow the given pattern at the given location. 
    void Ban(int location, int pattern)
    {
        // set that pattern at that location to false;
        outputArray[location][pattern] = false;

        // Zero out matching pattern count in each direction for this cell.
        int[] comp = matchingPatternsCount[location][pattern];
        for (int d = 0; d < 4; d++) comp[d] = 0;
        
        // put this location and pattern on the stack
        stack[stackPointer] = (location, pattern);
        stackPointer++;

        // update helper data
        sumsOfOnes[location] -= 1;
        sumsOfWeights[location] -= weights[pattern];
        sumsOfWeightLogWeights[location] -= weightLogWeights[pattern];

        double sum = sumsOfWeights[location];
        entropies[location] = Math.Log(sum) - sumsOfWeightLogWeights[location] / sum;
    }

    void Clear()
    {
        for (int i = 0; i < outputArray.Length; i++)
        {
            for (int t = 0; t < patternCount; t++)
            {
                outputArray[i][t] = true;
                for (int d = 0; d < 4; d++) matchingPatternsCount[i][t][d] = matchingPatternsInDir[opposite[d]][t].Length;
            }

            sumsOfOnes[i] = weights.Length;
            sumsOfWeights[i] = sumOfWeights;
            sumsOfWeightLogWeights[i] = sumOfWeightLogWeights;
            entropies[i] = startingEntropy;
            observed[i] = -1;
        }
        observedSoFar = 0;

        if (ground)
        {
            for (int x = 0; x < outputWidth; x++)
            {
                for (int t = 0; t < patternCount - 1; t++) Ban(x + (outputHeight - 1) * outputWidth, t);
                for (int y = 0; y < outputHeight - 1; y++) Ban(x + y * outputWidth, patternCount - 1);
            }
            Propagate();
        }
    }

    public abstract void Save(string filename);

    protected static int[] dx = { -1, 0, 1, 0 };
    protected static int[] dy = { 0, 1, 0, -1 };
    static int[] opposite = { 2, 3, 0, 1 };
}
