using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

public class MsaJob
{
    public string MsaName { get; set; }
    public string RefFasta { get; set; }
    public List<string> ExtractedFasta { get; set; }
}

public class MsaRunner
{
    private readonly string _referenceDir;
    private readonly string _extractedRootDir;
    private readonly string _msaRootDir;

    public MsaRunner(string referenceDir, string extractedRootDir, string msaRootDir)
    {
        _referenceDir = referenceDir;
        _extractedRootDir = extractedRootDir;
        _msaRootDir = msaRootDir;
        Directory.CreateDirectory(_msaRootDir);
    }

    public void RunAllForGenera(IEnumerable<string> genera)
    {
        var jobs = new List<MsaJob>
        {
            new MsaJob { MsaName = "mcya", RefFasta = "mcyA.fasta", ExtractedFasta = new List<string> { "mcyA.fasta", "adenylation.fasta" } },
            new MsaJob { MsaName = "mcyb", RefFasta = "mcyB.fasta", ExtractedFasta = new List<string> { "mcyB.fasta", "adenylation.fasta" } },
            new MsaJob { MsaName = "mcye", RefFasta = "mcyE.fasta", ExtractedFasta = new List<string> { "mcyE.fasta", "aminotransferase.fasta" } },
            new MsaJob { MsaName = "mcyh", RefFasta = "mcyH.fasta", ExtractedFasta = new List<string> { "mcyH.fasta", "ABC transporter.fasta" } }
        };

        foreach (var genus in genera)
        {
            string extractedDir = Path.Combine(_extractedRootDir, genus);
            string msaDir = Path.Combine(_msaRootDir, genus);
            Directory.CreateDirectory(msaDir);

            foreach (var job in jobs)
            {
                string msaInput = Path.Combine(msaDir, $"{job.MsaName}_msa.fasta");
                string msaOutput = Path.Combine(msaDir, $"{job.MsaName}_aligned.fasta");

                using (var writer = new StreamWriter(msaInput))
                {
                    string refPath = Path.Combine(_referenceDir, job.RefFasta);
                    if (File.Exists(refPath))
                    {
                        writer.WriteLine(File.ReadAllText(refPath).TrimEnd());
                    }
                    else
                    {
                        Console.WriteLine($"Warning: Reference file not found: {refPath}");
                    }

                    foreach (var extFasta in job.ExtractedFasta)
                    {
                        string extPath = Path.Combine(extractedDir, extFasta);
                        if (File.Exists(extPath))
                        {
                            writer.WriteLine(File.ReadAllText(extPath).TrimEnd());
                        }
                        else
                        {
                            Console.WriteLine($"Warning: Extracted file not found: {extPath}");
                        }
                    }
                }

                // Validate with external tool
                if (!ValidateFastaWithExternalTool(msaInput))
                {
                    Console.WriteLine($"Invalid FASTA format in {msaInput}. Skipping alignment.");
                    continue;
                }

                RunMafft(msaInput, msaOutput);
                Console.WriteLine($"[{genus}] MSA complete: {msaOutput}");
            }
        }
    }

    private void RunMafft(string inputFasta, string outputFasta)
    {
        var psi = new ProcessStartInfo
        {
            FileName = @"E:\Research\Research\muscle.exe",
            Arguments = $"-in \"{inputFasta}\" -out \"{outputFasta}\" -maxiters 1 -diags1",
            RedirectStandardOutput = false,
            RedirectStandardError = true,
            UseShellExecute = false,
            CreateNoWindow = true
        };
        using var process = Process.Start(psi);
        string error = process.StandardError.ReadToEnd();
        process.WaitForExit();
        if (!string.IsNullOrWhiteSpace(error))
            Console.WriteLine($"MUSCLE error: {error}");
    }

    private bool ValidateFastaWithExternalTool(string fastaPath)
    {
        string exePath = @"E:\Research\Research\Validate-Fasta-File\bin\ValidateFastaFile.exe";
        if (!File.Exists(exePath))
        {
            Console.WriteLine($"Error: Validation tool not found at {exePath}");
            return false;
        }

        var psi = new ProcessStartInfo
        {
            FileName = exePath,
            Arguments = $"\"{fastaPath}\"",
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            UseShellExecute = false,
            CreateNoWindow = true
        };

        try
        {
            using var process = Process.Start(psi);
            process.WaitForExit();

            string output = process.StandardOutput.ReadToEnd();
            string error = process.StandardError.ReadToEnd();

            if (!string.IsNullOrWhiteSpace(output))
                Console.WriteLine(output);
            if (!string.IsNullOrWhiteSpace(error))
                Console.WriteLine(error);

            return process.ExitCode == 0;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Failed to start validation tool: {ex.Message}");
            return false;
        }
    }
}