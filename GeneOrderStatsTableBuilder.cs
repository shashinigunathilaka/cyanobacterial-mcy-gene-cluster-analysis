using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

public class GeneOrderStats
{
    public string Gene { get; set; }
    public string Order { get; set; }
    public int PolymorphicSites { get; set; }
    public double NucleotideDiversity { get; set; }
    public double TajimasD { get; set; }
    public double RecombinationRate { get; set; }
    public double MutationRatePerBase { get; set; }
}

public static class GeneOrderStatsTableBuilder
{
    // List your genes and orders here
    private static readonly string[] Genes = { "mcyA", "mcyB", "mcyE", "mcyH" };
    private static readonly string[] Orders = { "Chroococcales", "Nostocales", "Oscillatoriales" };
    private static readonly string BasePath = @"E:\Research\Research\MSA";

    public static List<GeneOrderStats> BuildTable()
    {
        var statsList = new List<GeneOrderStats>();

        foreach (var gene in Genes)
        {
            foreach (var order in Orders)
            {
                string statsFile = Path.Combine(BasePath, order, $"{gene.ToLower()}_stats.tsv");
                if (!File.Exists(statsFile))
                    continue;

                var lines = File.ReadAllLines(statsFile);
                if (lines.Length < 2) continue;

                var header = lines[0].Split('\t');
                var values = lines[1].Split('\t');

                var stat = new GeneOrderStats
                {
                    Gene = gene,
                    Order = order,
                    PolymorphicSites = TryParseInt(GetValue(header, values, "PolymorphicSites")),
                    NucleotideDiversity = TryParseDouble(GetValue(header, values, "NucleotideDiversity")),
                    TajimasD = TryParseDouble(GetValue(header, values, "TajimasD")),
                    RecombinationRate = TryParseDouble(GetValue(header, values, "RecombinationRate")),
                    MutationRatePerBase = TryParseDouble(GetValue(header, values, "MutationRatePerBase"))
                };
                statsList.Add(stat);
            }
        }
        return statsList;
    }

    private static string GetValue(string[] header, string[] values, string colName)
    {
        int idx = Array.IndexOf(header, colName);
        return idx >= 0 ? values[idx] : "0";
    }

    private static int TryParseInt(string value)
    {
        return int.TryParse(value, out int result) ? result : 0;
    }

    private static double TryParseDouble(string value)
    {
        return double.TryParse(value, out double result) ? result : double.NaN;
    }

    public static void WriteTable(string outPath)
    {
        var stats = BuildTable();
        using var writer = new StreamWriter(outPath);
        writer.WriteLine("Gene\tOrder\tPolymorphicSites\tNucleotideDiversity\tTajimasD\tRecombinationRate\tMutationRatePerBase");
        foreach (var s in stats)
        {
            writer.WriteLine($"{s.Gene}\t{s.Order}\t{s.PolymorphicSites}\t{s.NucleotideDiversity}\t{s.TajimasD}\t{s.RecombinationRate}\t{s.MutationRatePerBase}");
        }
    }
}