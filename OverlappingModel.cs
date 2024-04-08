// Copyright (C) 2016 Maxim Gumin, The MIT License (MIT)

using System;
using System.Collections.Generic;

class OverlappingModel : Model
{
    List<byte[]> patterns;
    List<int> colors;

    public OverlappingModel(string name, int patternSize, int width, int height, bool periodicInput, bool periodic, int symmetry, bool ground, Heuristic heuristic)
        : base(width, height, patternSize, periodic, heuristic)
    {
        var (bitmap, SX, SY) = BitmapHelper.LoadBitmap($"samples/{name}.png");        
        byte[] sample = new byte[bitmap.Length];
        colors = new List<int>();
        for (int i = 0; i < sample.Length; i++)
        {
            int color = bitmap[i];
            int k = 0;
            for (; k < colors.Count; k++) if (colors[k] == color) break;
            if (k == colors.Count) colors.Add(color);
            sample[i] = (byte)k;
        }

        static byte[] pattern(Func<int, int, byte> f, int N)
        {
            byte[] result = new byte[N * N];
            for (int y = 0; y < N; y++) for (int x = 0; x < N; x++) result[x + y * N] = f(x, y);
            return result;
        };
        static byte[] rotate(byte[] p, int N) => pattern((x, y) => p[N - 1 - y + x * N], N);
        static byte[] reflect(byte[] p, int N) => pattern((x, y) => p[N - 1 - x + y * N], N);

        static long hash(byte[] p, int C)
        {
            long result = 0, power = 1;
            for (int i = 0; i < p.Length; i++)
            {
                result += p[p.Length - 1 - i] * power;
                power *= C;
            }
            return result;
        };

        patterns = new();
        Dictionary<long, int> patternIndices = new();
        List<double> weightList = new();

        int C = colors.Count;
        int xmax = periodicInput ? SX : SX - patternSize + 1;
        int ymax = periodicInput ? SY : SY - patternSize + 1;
        for (int y = 0; y < ymax; y++) for (int x = 0; x < xmax; x++)
            {
                byte[][] ps = new byte[8][];

                ps[0] = pattern((dx, dy) => sample[(x + dx) % SX + (y + dy) % SY * SX], patternSize);
                ps[1] = reflect(ps[0], patternSize);
                ps[2] = rotate(ps[0], patternSize);
                ps[3] = reflect(ps[2], patternSize);
                ps[4] = rotate(ps[2], patternSize);
                ps[5] = reflect(ps[4], patternSize);
                ps[6] = rotate(ps[4], patternSize);
                ps[7] = reflect(ps[6], patternSize);

                for (int k = 0; k < symmetry; k++)
                {
                    byte[] p = ps[k];
                    long h = hash(p, C);
                    if (patternIndices.TryGetValue(h, out int index)) weightList[index] = weightList[index] + 1;
                    else
                    {
                        patternIndices.Add(h, weightList.Count);
                        weightList.Add(1.0);
                        patterns.Add(p);
                    }
                }
            }

        weights = weightList.ToArray();
        patternCount = weights.Length;
        this.ground = ground;

        
        // p1 and p2: the two patterns to compare
        // dx/y = orthogonal offset
        // N = pattern size
        static bool agrees(byte[] p1, byte[] p2, int dx, int dy, int N)
        {
            int xmin = dx < 0 ? 0 : dx, xmax = dx < 0 ? dx + N : N, ymin = dy < 0 ? 0 : dy, ymax = dy < 0 ? dy + N : N;
            for (int y = ymin; y < ymax; y++) for (int x = xmin; x < xmax; x++) if (p1[x + N * y] != p2[x - dx + N * (y - dy)]) return false;
            return true;
        };

        // Propogator is a 3D array of ints
        // [Direction], [Pattern], [Pattern]
        // The array [D][P] is a list of the patterns that "agree" with P.
        // This is not a boolean vector; it's a ragged array of indices.
        matchingPatternsInDir = new int[4][][];
        for (int d = 0; d < 4; d++)
        {
            matchingPatternsInDir[d] = new int[patternCount][];
            for (int p1 = 0; p1 < patternCount; p1++)
            {
                List<int> list = new();
                for (int p2 = 0; p2 < patternCount; p2++)
                {
                    if (agrees(
                            patterns[p1],
                            patterns[p2],
                            dx[d],
                            dy[d],
                            patternSize))
                    {
                        list.Add(p2);
                    }
                }

                matchingPatternsInDir[d][p1] = new int[list.Count];
                for (int c = 0; c < list.Count; c++) matchingPatternsInDir[d][p1][c] = list[c];
            }
        }
    }

    public override void Save(string filename)
    {
        int[] bitmap = new int[outputWidth * outputHeight];
        if (observed[0] >= 0)
        {
            for (int y = 0; y < outputHeight; y++)
            {
                int dy = y < outputHeight - PatternSize + 1 ? 0 : PatternSize - 1;
                for (int x = 0; x < outputWidth; x++)
                {
                    int dx = x < outputWidth - PatternSize + 1 ? 0 : PatternSize - 1;
                    bitmap[x + y * outputWidth] = colors[patterns[observed[x - dx + (y - dy) * outputWidth]][dx + dy * PatternSize]];
                }
            }
        }
        else
        {
            for (int i = 0; i < outputArray.Length; i++)
            {
                int contributors = 0, r = 0, g = 0, b = 0;
                int x = i % outputWidth, y = i / outputWidth;
                for (int dy = 0; dy < PatternSize; dy++) for (int dx = 0; dx < PatternSize; dx++)
                    {
                        int sx = x - dx;
                        if (sx < 0) sx += outputWidth;

                        int sy = y - dy;
                        if (sy < 0) sy += outputHeight;

                        int s = sx + sy * outputWidth;
                        if (!periodic && (sx + PatternSize > outputWidth || sy + PatternSize > outputHeight || sx < 0 || sy < 0)) continue;
                        for (int t = 0; t < patternCount; t++) if (outputArray[s][t])
                            {
                                contributors++;
                                int argb = colors[patterns[t][dx + dy * PatternSize]];
                                r += (argb & 0xff0000) >> 16;
                                g += (argb & 0xff00) >> 8;
                                b += argb & 0xff;
                            }
                    }
                bitmap[i] = unchecked((int)0xff000000 | ((r / contributors) << 16) | ((g / contributors) << 8) | b / contributors);
            }
        }
        BitmapHelper.SaveBitmap(bitmap, outputWidth, outputHeight, filename);
    }
}
