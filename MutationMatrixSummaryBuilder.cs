using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

public static class MutationMatrixSummaryBuilder
{
    public static string MatrixDir = @"E:\Research\Research\MSA\feature_matrices";
    public static string OutputFile = @"E:\Research\Research\MSA\feature_matrices\mutation_matrix_summary.csv";

    public static void BuildSummary()
    {
        var summaryRows = new List<string>();
        summaryRows.Add("Gene,Order,SeqHeader,PercentMutated,NumMutated,NumPositions");

        var files = Directory.GetFiles(MatrixDir, "*_mutation_matrix.csv");
        foreach (var file in files)
        {
            var fname = Path.GetFileName(file);
            var parts = fname.Split('_');
            if (parts.Length < 3) continue;
            string order = parts[0];
            string gene = parts[1];

            using var reader = new StreamReader(file);
            var headerLine = reader.ReadLine();
            if (headerLine == null) continue;
            var colHeaders = headerLine.Split(',').Skip(1).ToArray(); // skip SeqHeader
            int numPositions = colHeaders.Length;

            string line;
            while ((line = reader.ReadLine()) != null)
            {
                var tokens = line.Split(',');
                if (tokens.Length < numPositions + 1) continue;
                string seqHeader = tokens[0];
                int numMutated = tokens.Skip(1).Count(x => x == "1");
                double percent = 100.0 * numMutated / numPositions;
                summaryRows.Add($"{gene},{order},{seqHeader},{percent:F2},{numMutated},{numPositions}");
            }
        }

        File.WriteAllLines(OutputFile, summaryRows);
        Console.WriteLine($"Summary written to: {OutputFile}");
    }
}