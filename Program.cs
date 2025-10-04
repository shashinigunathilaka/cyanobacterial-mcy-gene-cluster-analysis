using System;
using System.IO;
using System.Text.RegularExpressions;
using System.Collections.Generic;

class Program
{
    // Genes of interest
    static readonly string[] TargetGenes = { "mcyA", "mcyB", "mcyE", "mcyH", "adenylation", "ABC transporter", "aminotransferase" };

    // Output structure: Genus -> Gene -> List of (Header, Sequence)
    static readonly Dictionary<string, Dictionary<string, List<(string, string)>>> GeneSequences =
        new(StringComparer.OrdinalIgnoreCase);

    static void Main()
    {
        // Set absolute paths for input and output
        string root = @"E:\Research\Research\DataSet";
        string referenceRoot = @"E:\Research\Research\ReferenceGenes";
        string outputRoot = @"E:\Research\Research\ExtractedGenes";

        if (!Directory.Exists(root))
        {
            Console.WriteLine($"DataSet directory not found at: {root}");
            return;
        }

        // -------------------------------------------------------------------------
        // The following block is commented out to prevent re-extraction of genes
        // and re-generation of FASTA files, which can be time-consuming and is
        // unnecessary if you already have the required files from a previous run.
        // Uncomment if you need to regenerate the extracted FASTA files.
        
        // Process all genus directories
        foreach (var genusDir in Directory.GetDirectories(root))
        {
            string genus = Path.GetFileName(genusDir);
            TraverseGenus(genus, genusDir);
        }

        // Process ReferenceGenes as a special genus
        if (Directory.Exists(referenceRoot))
        {
            TraverseGenus("ReferenceGenes", referenceRoot);
        }

        // Output summary and save results
        foreach (var genus in GeneSequences.Keys)
        {
            foreach (var gene in GeneSequences[genus].Keys)
            {
                string outDir = Path.Combine(outputRoot, genus);
                Directory.CreateDirectory(outDir);
                string outFile = Path.Combine(outDir, $"{gene}.fasta");
                using var writer = new StreamWriter(outFile);
                WriteFastaSequencesUnique(writer, GeneSequences[genus][gene]);
                Console.WriteLine($"[{genus}] {gene}: {GeneSequences[genus][gene].Count} sequences saved to {outFile}");
            }
        }
        
        // -------------------------------------------------------------------------

        // -------------------------------------------------------------------------
        // The following block is commented out to prevent re-running the alignment
        // process (MSA), which is computationally expensive and unnecessary if
        // alignments have already been generated. Uncomment if you need to realign.
        
        var genera = new[] { "Chroococcales", "Nostocales", "Oscillatoriales" };
        var msaRunner = new MsaRunner(
            @"E:\Research\Research\ReferenceGenes", // Correct path
            outputRoot,
            @"E:\Research\Research\MSA"
        );
        msaRunner.RunAllForGenera(genera);
        
        // -------------------------------------------------------------------------

        // --- Mutation Analysis on existing alignment files ---
        // This section will run mutation analysis on the alignment files that were
        // previously generated, without repeating extraction or alignment steps.
        
        var jobs = new[] { "mcya", "mcyb", "mcye", "mcyh" }; // Use the jobs you actually aligned

        foreach (var genus in genera)
        {
            foreach (var job in jobs)
            {
                string alignedFastaPath = $@"E:\Research\Research\MSA\{genus}\{job}_aligned.fasta";
                if (!File.Exists(alignedFastaPath))
                {
                    Console.WriteLine($"Alignment file not found: {alignedFastaPath}");
                    continue;
                }
                string referenceHeader = GetReferenceHeader(alignedFastaPath);
                string summaryPath = $@"E:\Research\Research\MSA\{genus}\{job}_mutation_summary.tsv";
                MutationAnalyzer.WriteMutationSummaryByOrder(alignedFastaPath, referenceHeader, summaryPath);
            }
        }


        //FeatureMatrixOneHotBuilder.BuildAllOneHotMatrices();
        FeatureMatrixBuilder.BuildMutationTypeSummary();
        MutationMatrixSummaryBuilder.BuildSummary();
        GeneOrderStatsTableBuilder.WriteTable(@"E:\Research\Research\MSA\plots\gene_order_stats_table.tsv");
        Console.WriteLine("Gene-order stats table written to plots folder.");
    }

    static void TraverseGenus(string genus, string genusDir)
    {
        foreach (var subDir in Directory.GetDirectories(genusDir, "*", SearchOption.AllDirectories))
        {
            foreach (var cdsFile in Directory.GetFiles(subDir, "cds_from_genomic.fna"))
            {
                ExtractGenesFromFasta(genus, cdsFile);
            }
        }
    }

    static void ExtractGenesFromFasta(string genus, string fastaPath)
    {
        string[] lines = File.ReadAllLines(fastaPath);
        string header = null;
        var seqBuilder = new System.Text.StringBuilder();
        string foundGene = null;

        foreach (var line in lines)
        {
            if (line.StartsWith(">"))
            {
                // Save previous sequence if it matches
                if (header != null && foundGene != null)
                {
                    AddGeneSequence(genus, foundGene, header, seqBuilder.ToString());
                }
                header = line.Substring(1).Trim();
                seqBuilder.Clear();
                foundGene = null;

                // Try to find gene name in header
                foreach (var gene in TargetGenes)
                {
                    // Look for [gene=...] or [protein=...] containing the gene keyword
                    if (Regex.IsMatch(header, $@"\[gene\s*=\s*{Regex.Escape(gene)}\]", RegexOptions.IgnoreCase) ||
                        Regex.IsMatch(header, $@"\[protein\s*=\s*[^\]]*{Regex.Escape(gene)}[^\]]*\]", RegexOptions.IgnoreCase))
                    {
                        foundGene = gene;
                        break;
                    }
                }
            }
            else if (header != null)
            {
                seqBuilder.Append(line.Trim());
            }
        }
        // Save last sequence
        if (header != null && foundGene != null)
        {
            AddGeneSequence(genus, foundGene, header, seqBuilder.ToString());
        }
    }

    static void AddGeneSequence(string genus, string gene, string header, string sequence)
    {
        if (!GeneSequences.ContainsKey(genus))
            GeneSequences[genus] = new Dictionary<string, List<(string, string)>>(StringComparer.OrdinalIgnoreCase);
        if (!GeneSequences[genus].ContainsKey(gene))
            GeneSequences[genus][gene] = new List<(string, string)>();
        GeneSequences[genus][gene].Add((header, sequence));
    }

    private static string ExtractProteinName(string header)
    {
        // Example: lcl|CP073041.1_cds_UXE61484.1_33 [protein=SomeProteinName]
        var match = Regex.Match(header, @"\[protein\s*=\s*([^\]]+)\]", RegexOptions.IgnoreCase);
        if (match.Success)
            return match.Groups[1].Value.Trim();
        // Fallback: use the whole header if no protein name found
        return header;
    }

    static void WriteFastaSequencesUnique(StreamWriter writer, List<(string header, string seq)> entries, int maxSequences = 1000)
    {
        var seen = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        int count = 0;
        foreach (var (header, seq) in entries)
        {
            if (count >= maxSequences) break;
            var proteinName = ExtractProteinName(header);
            if (!seen.Add(proteinName))
            {
                Console.WriteLine($"Warning: Duplicate protein name found: {proteinName}, skipping.");
                continue;
            }
            WriteFastaSequence(writer, header, seq);
            count++;
        }
    }

    static void WriteFastaSequence(StreamWriter writer, string header, string sequence, int lineLength = 60)
    {
        writer.WriteLine($">{header}");
        for (int i = 0; i < sequence.Length; i += lineLength)
        {
            int len = Math.Min(lineLength, sequence.Length - i);
            writer.WriteLine(sequence.Substring(i, len));
        }
    }

    // Helper to extract the first header as reference (customize as needed)
    static string GetReferenceHeader(string fastaPath)
    {
        foreach (var line in File.ReadLines(fastaPath))
        {
            if (line.StartsWith(">"))
                return line.Substring(1).Trim();
        }
        return null;
    }
}
