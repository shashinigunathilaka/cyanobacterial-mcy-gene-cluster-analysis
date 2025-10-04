using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

public static class FeatureMatrixOneHotBuilder
{
    public static readonly string[] Genes = { "mcyA", "mcyB", "mcyE", "mcyH" };
    public static readonly string[] Orders = { "Chroococcales", "Nostocales", "Oscillatoriales" };
    public static readonly char[] Bases = { 'A', 'C', 'G', 'T', '-' };

    public static string BasePath = @"E:\Research\Research\MSA";
    public static string OutputDir = Path.Combine(BasePath, "feature_matrices_onehot");

    public static void BuildAllOneHotMatrices()
    {
        Directory.CreateDirectory(OutputDir);

        foreach (var gene in Genes)
        {
            foreach (var order in Orders)
            {
                string fastaPath = Path.Combine(BasePath, order, $"{gene.ToLower()}_aligned.fasta");
                if (!File.Exists(fastaPath))
                {
                    Console.WriteLine($"File not found: {fastaPath}");
                    continue;
                }

                var matrix = BuildOneHotMatrix(fastaPath, out var headers, out var colNames);
                if (matrix == null) continue;

                string outFile = Path.Combine(OutputDir, $"{order}_{gene}_onehot_matrix.csv");
                WriteMatrixToCsv(matrix, headers, colNames, outFile);
                Console.WriteLine($"Saved: {outFile}");
            }
        }
    }

    // Returns: rows = sequence headers, columns = positions*bases, values = 0/1
    public static int[,] BuildOneHotMatrix(string fastaPath, out List<string> headers, out List<string> colNames)
    {
        var seqs = MutationAnalyzer.ParseAlignedFasta(fastaPath);
        headers = seqs.Keys.ToList();
        if (headers.Count == 0) { colNames = null; return null; }

        int seqCount = headers.Count;
        int seqLen = seqs[headers[0]].Length;
        int baseCount = Bases.Length;
        int[,] matrix = new int[seqCount, seqLen * baseCount];

        colNames = new List<string>();
        for (int pos = 0; pos < seqLen; pos++)
        {
            foreach (var b in Bases)
                colNames.Add($"Pos{pos + 1}_{b}");
        }

        for (int i = 0; i < seqCount; i++)
        {
            var seq = seqs[headers[i]];
            for (int pos = 0; pos < seqLen; pos++)
            {
                char baseAtPos = (pos < seq.Length) ? seq[pos] : '-';
                for (int b = 0; b < baseCount; b++)
                {
                    matrix[i, pos * baseCount + b] = (baseAtPos == Bases[b]) ? 1 : 0;
                }
            }
        }
        return matrix;
    }

    public static void WriteMatrixToCsv(int[,] matrix, List<string> rowHeaders, List<string> colHeaders, string outFile)
    {
        using var writer = new StreamWriter(outFile);
        // Write the first column header explicitly
        writer.WriteLine("SeqHeader," + string.Join(",", colHeaders));
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            writer.Write(rowHeaders[i]);
            for (int j = 0; j < matrix.GetLength(1); j++)
                writer.Write($",{matrix[i, j]}");
            writer.WriteLine();
        }
    }
}