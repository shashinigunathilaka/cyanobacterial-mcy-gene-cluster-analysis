using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

public static class FeatureMatrixBuilder
{
    // List of genes and orders to process
    public static readonly string[] Genes = { "mcyA", "mcyB", "mcyE", "mcyH" };
    public static readonly string[] Orders = { "Chroococcales", "Nostocales", "Oscillatoriales" };

    // Base path for aligned FASTA files and output
    public static string BasePath = @"E:\Research\Research\MSA";
    public static string OutputDir = Path.Combine(BasePath, "feature_matrices");

    public static void BuildAllFeatureMatrices()
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

                var matrix = BuildMutationMatrix(fastaPath);
                if (matrix == null) continue;

                string outFile = Path.Combine(OutputDir, $"{order}_{gene}_mutation_matrix.csv");
                WriteMatrixToCsv(matrix, outFile);
                Console.WriteLine($"Saved: {outFile}");
            }
        }
    }

    // Returns: rows = sequence headers, columns = positions, values = 0/1
    public static int[,] BuildMutationMatrix(string fastaPath)
    {
        var seqs = MutationAnalyzer.ParseAlignedFasta(fastaPath);
        if (seqs.Count == 0) return null;

        var headers = seqs.Keys.ToList();
        var reference = seqs[headers[0]];
        int seqCount = headers.Count;
        int seqLen = reference.Length;
        int[,] matrix = new int[seqCount, seqLen];

        for (int i = 0; i < seqCount; i++)
        {
            var seq = seqs[headers[i]];
            for (int j = 0; j < seqLen; j++)
            {
                if (j >= seq.Length || j >= reference.Length)
                {
                    matrix[i, j] = 1; // treat as mutation if missing
                }
                else
                {
                    matrix[i, j] = (seq[j] == reference[j]) ? 0 : 1;
                }
            }
        }
        return matrix;
    }

    // Writes the matrix to CSV with headers
    public static void WriteMatrixToCsv(int[,] matrix, string outFile)
    {
        using var writer = new StreamWriter(outFile);
        int rowCount = matrix.GetLength(0);
        int colCount = matrix.GetLength(1);

        // Write column headers
        var colHeaders = Enumerable.Range(1, colCount).Select(i => $"Pos{i}");
        writer.WriteLine("SeqHeader," + string.Join(",", colHeaders));

        // Write each row
        for (int i = 0; i < rowCount; i++)
        {
            writer.Write($"Seq{i + 1}");
            for (int j = 0; j < colCount; j++)
            {
                writer.Write($",{matrix[i, j]}");
            }
            writer.WriteLine();
        }
    }

    public static void BuildMutationTypeSummary()
    {
        var summaryRows = new List<string>();
        summaryRows.Add("Gene,Order,Missense,Frameshift,FrameshiftDeletion,FrameshiftInsertion,Complex,InFrame");

        foreach (var gene in Genes)
        {
            foreach (var order in Orders)
            {
                string fastaPath = Path.Combine(BasePath, order, $"{gene.ToLower()}_aligned.fasta");
                if (!File.Exists(fastaPath)) continue;

                var seqs = MutationAnalyzer.ParseAlignedFasta(fastaPath);
                if (seqs.Count == 0) continue;

                var headers = seqs.Keys.ToList();
                var reference = seqs[headers[0]];
                // Initialize counters
                var typeCounts = new Dictionary<MutationAnalyzer.MutationType, int>();
                foreach (MutationAnalyzer.MutationType mt in Enum.GetValues(typeof(MutationAnalyzer.MutationType)))
                    typeCounts[mt] = 0;

                // For each sequence (except reference)
                for (int i = 1; i < headers.Count; i++)
                {
                    var seq = seqs[headers[i]];
                    for (int j = 0; j < reference.Length; j++)
                    {
                        if (j >= seq.Length) continue;
                        var mtype = MutationAnalyzer.ClassifyMutation(reference, seq, j);
                        if (mtype != MutationAnalyzer.MutationType.None)
                            typeCounts[mtype]++;
                    }
                }

                summaryRows.Add($"{gene},{order},{typeCounts[MutationAnalyzer.MutationType.Missense]},{typeCounts[MutationAnalyzer.MutationType.Frameshift]},{typeCounts[MutationAnalyzer.MutationType.FrameshiftDeletion]},{typeCounts[MutationAnalyzer.MutationType.FrameshiftInsertion]},{typeCounts[MutationAnalyzer.MutationType.Complex]},{typeCounts[MutationAnalyzer.MutationType.InFrame]}");
            }
        }

        string outFile = Path.Combine(OutputDir, "mutation_type_summary.csv");
        File.WriteAllLines(outFile, summaryRows);
        Console.WriteLine($"Mutation type summary written to: {outFile}");
    }
}