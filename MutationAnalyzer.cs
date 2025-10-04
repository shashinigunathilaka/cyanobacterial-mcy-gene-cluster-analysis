using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

public class MutationAnalyzer
{
    public enum MutationType
    {
        None,
        Missense,
        Frameshift,
        FrameshiftDeletion,
        FrameshiftInsertion,
        Complex,
        InFrame
    }

    // Parse aligned FASTA into a dictionary: header -> sequence
    public static Dictionary<string, string> ParseAlignedFasta(string path)
    {
        var dict = new Dictionary<string, string>();
        string currentHeader = null;
        var seqBuilder = new StringBuilder();

        foreach (var line in File.ReadLines(path))
        {
            if (line.StartsWith(">"))
            {
                if (currentHeader != null)
                    dict[currentHeader] = seqBuilder.ToString();
                currentHeader = line.Substring(1).Trim();
                seqBuilder.Clear();
            }
            else if (currentHeader != null)
            {
                seqBuilder.Append(line.Trim());
            }
        }
        if (currentHeader != null)
            dict[currentHeader] = seqBuilder.ToString();

        return dict;
    }

    // Count mutations (mismatches, ignoring gaps) between reference and another sequence
    public static int CountMutations(string reference, string sequence)
    {
        int mutations = 0;
        for (int i = 0; i < Math.Min(reference.Length, sequence.Length); i++)
        {
            char refChar = reference[i];
            char seqChar = sequence[i];
            if (refChar == '-' || seqChar == '-') continue; // ignore gaps
            if (refChar != seqChar) mutations++;
        }
        return mutations;
    }

    // Example: extract genus from header (customize as needed)
    public static string ExtractGenus(string header)
    {
        // Example: genus is first word in header, or use regex if needed
        var match = Regex.Match(header, @"\b([A-Za-z]+)\b");
        return match.Success ? match.Groups[1].Value : "Unknown";
    }

    // Extract order from header (example: second word in header)
    public static string ExtractOrder(string header)
    {
        var match = Regex.Match(header, @"\[order\s*=\s*([^\]]+)\]", RegexOptions.IgnoreCase);
        if (match.Success)
            return match.Groups[1].Value.Trim();
        return "Unknown";
    }

    private static string ExtractOrderFromPath(string fastaPath)
    {
        var dir = Path.GetDirectoryName(fastaPath);
        if (string.IsNullOrEmpty(dir)) return "Unknown";
        var order = Path.GetFileName(dir);
        return string.IsNullOrEmpty(order) ? "Unknown" : order;
    }

    // Summarize mutation counts by genus
    public static void SummarizeMutations(string alignedFastaPath, string referenceHeader)
    {
        var seqs = ParseAlignedFasta(alignedFastaPath);
        if (!seqs.TryGetValue(referenceHeader, out var reference))
        {
            Console.WriteLine($"Reference header '{referenceHeader}' not found in alignment.");
            return;
        }

        var genusCounts = new Dictionary<string, List<int>>();

        foreach (var kvp in seqs)
        {
            if (kvp.Key == referenceHeader) continue;
            int muts = CountMutations(reference, kvp.Value);
            string genus = ExtractGenus(kvp.Key);
            if (!genusCounts.ContainsKey(genus))
                genusCounts[genus] = new List<int>();
            genusCounts[genus].Add(muts);
            Console.WriteLine($"Seq: {kvp.Key} | Genus: {genus} | Mutations: {muts}");
        }

        Console.WriteLine("\nSummary by genus:");
        foreach (var kvp in genusCounts)
        {
            double avg = kvp.Value.Count > 0 ? kvp.Value.Average() : 0;
            Console.WriteLine($"{kvp.Key}: {kvp.Value.Count} sequences, Avg mutations: {avg:F2}");
        }
    }

    // Summarize mutation counts by order
    public static void SummarizeMutationsByOrder(string alignedFastaPath, string referenceHeader)
    {
        var seqs = ParseAlignedFasta(alignedFastaPath);
        if (!seqs.TryGetValue(referenceHeader, out var reference))
        {
            Console.WriteLine($"Reference header '{referenceHeader}' not found in alignment.");
            return;
        }

        var orderCounts = new Dictionary<string, List<int>>();

        foreach (var kvp in seqs)
        {
            if (kvp.Key == referenceHeader) continue;
            int muts = CountMutations(reference, kvp.Value);
            string order = ExtractOrder(kvp.Key);
            if (!orderCounts.ContainsKey(order))
                orderCounts[order] = new List<int>();
            orderCounts[order].Add(muts);
            Console.WriteLine($"Seq: {kvp.Key} | Order: {order} | Mutations: {muts}");
        }

        Console.WriteLine("\nSummary by order:");
        foreach (var kvp in orderCounts)
        {
            double avg = kvp.Value.Count > 0 ? kvp.Value.Average() : 0;
            Console.WriteLine($"{kvp.Key}: {kvp.Value.Count} sequences, Avg mutations: {avg:F2}");
        }
    }

    // Write mutation summary to a file
    public static void WriteMutationSummary(
        string alignedFastaPath,
        string referenceHeader,
        string outputSummaryPath)
    {
        var seqs = ParseAlignedFasta(alignedFastaPath);
        if (!seqs.TryGetValue(referenceHeader, out var reference))
        {
            Console.WriteLine($"Reference header '{referenceHeader}' not found in alignment.");
            return;
        }

        var genusCounts = new Dictionary<string, List<int>>();
        var perSequence = new List<(string Header, string Genus, int Mutations)>();

        foreach (var kvp in seqs)
        {
            if (kvp.Key == referenceHeader) continue;
            int muts = CountMutations(reference, kvp.Value);
            string genus = ExtractGenus(kvp.Key);
            if (!genusCounts.ContainsKey(genus))
                genusCounts[genus] = new List<int>();
            genusCounts[genus].Add(muts);
            perSequence.Add((kvp.Key, genus, muts));
        }

        using (var writer = new StreamWriter(outputSummaryPath))
        {
            writer.WriteLine("Header\tGenus\tMutations");
            foreach (var entry in perSequence)
                writer.WriteLine($"{entry.Header}\t{entry.Genus}\t{entry.Mutations}");

            writer.WriteLine();
            writer.WriteLine("Genus\tNumSequences\tAvgMutations");
            foreach (var kvp in genusCounts)
            {
                double avg = kvp.Value.Count > 0 ? kvp.Value.Average() : 0;
                writer.WriteLine($"{kvp.Key}\t{kvp.Value.Count}\t{avg:F2}");
            }
        }
        Console.WriteLine($"Mutation summary written to {outputSummaryPath}");
    }

    // Write mutation summary to a file (by order)
    public static void WriteMutationSummaryByOrder(
        string alignedFastaPath,
        string referenceHeader,
        string outputSummaryPath)
    {
        var seqs = ParseAlignedFasta(alignedFastaPath);
        if (!seqs.TryGetValue(referenceHeader, out var reference))
        {
            Console.WriteLine($"Reference header '{referenceHeader}' not found in alignment.");
            return;
        }

        // Use folder name as order
        string order = ExtractOrderFromPath(alignedFastaPath);

        var orderCounts = new Dictionary<string, List<int>>();
        var perSequence = new List<(string Header, string Order, int Mutations)>();

        foreach (var kvp in seqs)
        {
            if (kvp.Key == referenceHeader) continue;
            int muts = CountMutations(reference, kvp.Value);
            // Use the extracted order for all sequences in this file
            if (!orderCounts.ContainsKey(order))
                orderCounts[order] = new List<int>();
            orderCounts[order].Add(muts);
            perSequence.Add((kvp.Key, order, muts));
        }

        using (var writer = new StreamWriter(outputSummaryPath))
        {
            writer.WriteLine("Header\tOrder\tMutations");
            foreach (var entry in perSequence)
                writer.WriteLine($"{entry.Header}\t{entry.Order}\t{entry.Mutations}");

            writer.WriteLine();
            writer.WriteLine("Order\tNumSequences\tAvgMutations");
            foreach (var kvp in orderCounts)
            {
                double avg = kvp.Value.Count > 0 ? kvp.Value.Average() : 0;
                writer.WriteLine($"{kvp.Key}\t{kvp.Value.Count}\t{avg:F2}");
            }
        }
        Console.WriteLine($"Mutation summary by order written to {outputSummaryPath}");
    }

    // Example: classify mutation type at a position (stub logic)
    public static MutationType ClassifyMutation(string reference, string sequence, int pos)
    {
        char refChar = reference[pos];
        char seqChar = sequence[pos];
        if (refChar == seqChar) return MutationType.None;
        if (seqChar == '-') return MutationType.FrameshiftDeletion;
        if (refChar == '-') return MutationType.FrameshiftInsertion;
        return MutationType.Missense;
    }
}