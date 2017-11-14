using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ASG2
{
    class Program
    {
        public static DateTime GetLinkerTime(TimeZoneInfo target = null)
        {
            var filePath = System.Reflection.Assembly.GetExecutingAssembly().Location;
            const int c_PeHeaderOffset = 60;
            const int c_LinkerTimestampOffset = 8;

            var buffer = new byte[2048];

            using (var stream = new FileStream(filePath, FileMode.Open, FileAccess.Read))
                stream.Read(buffer, 0, 2048);

            var offset = BitConverter.ToInt32(buffer, c_PeHeaderOffset);
            var secondsSince1970 = BitConverter.ToInt32(buffer, offset + c_LinkerTimestampOffset);
            var epoch = new DateTime(1970, 1, 1, 0, 0, 0, DateTimeKind.Utc);

            var linkTimeUtc = epoch.AddSeconds(secondsSince1970);

            var tz = target ?? TimeZoneInfo.Local;
            var localTime = TimeZoneInfo.ConvertTimeFromUtc(linkTimeUtc, tz);

            return localTime;
        }

        static int min_peptide_length = 6;
        static int max_peptide_length = 50;
        static float max_peptide_mass = 4600;

        static List<string> skipList = new List<string>();

        static int peptidePerGeneLimit;

        static FastaParser fastaParser;

        static void Main(string[] args)
        {
            Console.Error.WriteLine("build time: " + GetLinkerTime());

            if (args.Length != 6)
            {
                Console.WriteLine("args length error!");
                Console.ReadLine();
                return;
            }

            int limit_unused = int.Parse(args[4]);
            peptidePerGeneLimit = int.Parse(args[5]);

            string fastaPath = Path.GetFullPath(args[0]);
            string gffPath = Path.GetFullPath(args[1]);
            string prefix = args[2];
            string outputPath = Path.GetFullPath(args[3]);

            fastaParser = new FastaParser(fastaPath);
            GffParser gffParser = new GffParser(new string[] { "biotype" });
            List<GffParser.Node> gffNodes = gffParser.ReadAllNodes(gffPath);
            gffNodes.RemoveAll((node) => node.feature == "chromosome" || node.feature == "supercontig" || node.seqid == "MT");
            foreach (var gene in gffNodes)
            {
                gene.subNodes.RemoveAll((transcript) => transcript.attributes["biotype"] != "protein_coding" && transcript.attributes["biotype"] != "polymorphic_pseudogene");
            }
            gffNodes.RemoveAll((node) => node.subNodes.Count == 0);

            List<Task> taskList = new List<Task>();
            foreach (var gene in gffNodes)
            {
                Task<string> readTask = Task.Run(() => fastaParser.ReadBlock(gene.seqid, gene.start, gene.end, gene.strand == Strand.Minus));
                List<int> startList = new List<int>();
                List<int> donorList = new List<int>();
                List<int> acceptorList = new List<int>();
                List<int> endingList = new List<int>();
                foreach (var transcript in gene.subNodes)
                {
                    int start;
                    int[] donor;
                    int[] acceptor;
                    int ending;
                    GetTranscriptInfo(transcript, out start, out donor, out acceptor, out ending);
                    startList.Add(start);
                    donorList.AddRange(donor);
                    acceptorList.AddRange(acceptor);
                    endingList.Add(ending);
                }
                int[] startArray = startList.Distinct().ToArray();
                int[] donorArray = donorList.Distinct().ToArray();
                int[] acceptorArray = acceptorList.Distinct().ToArray();
                int[] endingArray = endingList.Distinct().ToArray();
                int offset = gene.start - 1;
                if (gene.strand == Strand.Minus)
                {
                    int length = fastaParser.GetSequenceLength(gene.seqid);
                    ReverseLocation(startArray, length);
                    ReverseLocation(donorArray, length);
                    ReverseLocation(acceptorArray, length);
                    ReverseLocation(endingArray, length);
                    offset = length - gene.end;
                }
                readTask.Wait();
                Sequence sequence = new Sequence(readTask.Result, offset);
                string reference = gene.ID.Split(':')[1];
                if (limit_unused != -1 && acceptorArray.Length > limit_unused)
                {
                    skipList.Add(reference);
                    continue;
                }
                //Work(reference, sequence, startArray, donorArray, acceptorArray, endingArray);
                Task task = Task.Run(() => Work(gene.seqid, gene.strand, reference, sequence, startArray, donorArray, acceptorArray, endingArray));
                taskList.Add(task);
            }
            gffNodes = null;
            Task[] taskArray = taskList.ToArray();
            //while (true)
            //{
            //    Task.WaitAny(taskArray, 1000);
            //    int count = taskArray.Where((t) => !t.IsCompleted).Count();
            //    Console.WriteLine(count);
            //    if (count == 0)
            //    {
            //        break;
            //    }
            //}
            Task.WaitAll(taskArray);
            foreach (var task in taskArray)
            {
                if (task.IsCompleted)
                {
                    var result = (task as Task<Tuple<string, string, Dictionary<string, HashSet<int[]>>>>).Result;
                    string reference = result.Item2;
                    Dictionary<string, HashSet<int[]>> peptideList = result.Item3;
                    foreach (var peptideKV in peptideList)
                    {
                        string peptide = peptideKV.Key;
                        if (!peptideSet.ContainsKey(peptide))
                        {
                            peptideSet.Add(peptide, new List<PeptideInfo>());
                        }
                        List<PeptideInfo> list = peptideSet[peptide];
                        foreach (var region in peptideKV.Value)
                        {
                            PeptideInfo info = new PeptideInfo { seq_id = result.Item1, reference = reference, region = region };
                            peptideSet[peptide].Add(info);
                        }
                    }
                    taskList.Remove(task);
                }
            }

            OutputResult(outputPath, prefix);
        }

        static bool CheckSpliceSite(Sequence sequence, int location, int type)
        {
            if (type == 1)
            {
                if (location + 1 >= sequence.Start && location + 2 <= sequence.End && sequence.GetSequence(location + 1, 2) == "GT")
                {
                    return true;
                }
                if (location + 1 >= sequence.Start && location + 2 <= sequence.End && sequence.GetSequence(location + 1, 2) == "GC")
                {
                    return true;
                }
            }
            if (type == 2)
            {
                if (location - 2 >= sequence.Start && location - 1 <= sequence.End && sequence.GetSequence(location - 2, 2) == "AG")
                {
                    return true;
                }
            }
            return false;
        }

        static bool GetTranscriptInfo(GffParser.Node transcript, out int start, out int[] donor, out int[] acceptor, out int ending)
        {
            var cds = transcript.subNodes.Where((node) => node.feature == "CDS");
            bool found = false;
            byte phase = 0;
            if (transcript.strand == Strand.Minus)
            {
                start = 0;
                ending = int.MaxValue;
                foreach (var node in cds)
                {
                    if (node.end > start)
                    {
                        if (node.end - phase >= node.start)
                        {
                            start = node.end;
                            phase = node.phase.Value;
                            found = true;
                        }
                    }
                    if (node.start < ending)
                    {
                        ending = node.start;
                    }
                }
                start -= phase;
            }
            else
            {
                start = int.MaxValue;
                ending = 0;
                foreach (var node in cds)
                {
                    if (node.start < start)
                    {
                        if (node.start + phase <= node.end)
                        {
                            start = node.start;
                            phase = node.phase.Value;
                            found = true;
                        }
                    }
                    if (node.end > ending)
                    {
                        ending = node.end;
                    }
                }
                start += phase;
            }

            if (!found)
            {
                throw new Exception();
            }

            List<int> location = new List<int>();
            foreach (var node in transcript.subNodes)
            {
                if (node.feature == "exon")
                {
                    location.Add(node.start);
                    location.Add(node.end);
                }
            }
            if (transcript.strand == Strand.Minus)
            {
                location.Sort((x, y) => y - x);
            }
            else
            {
                location.Sort((x, y) => x - y);
            }
            List<int> donorList = new List<int>();
            List<int> acceptorList = new List<int>();
            for (int i = 1; i < location.Count - 1; i++)
            {
                if (i % 2 == 0)
                {
                    acceptorList.Add(location[i]);
                }
                else
                {
                    donorList.Add(location[i]);
                }
            }
            donor = donorList.ToArray();
            acceptor = acceptorList.ToArray();
            return true;
        }

        static int[] ReverseLocation(int[] source, int length)
        {
            for (int i = 0; i < source.Length; i++)
            {
                source[i] = length - source[i] + 1;
            }
            return source;
        }

        struct PeptideInfo : IEquatable<PeptideInfo>
        {
            public string seq_id;
            public string reference;
            public int[] region;

            public bool Equals(PeptideInfo other)
            {
                if (seq_id == other.seq_id && reference == other.reference)
                {
                    return region.SequenceEqual(other.region);
                }
                return false;
            }
        }

        static Dictionary<string, List<PeptideInfo>> peptideSet = new Dictionary<string, List<PeptideInfo>>();

        // 未使用
        static void AddPeptide(string peptide, string reference, int[] region)
        {
            if (peptide.Length < 5 || peptide.Length > 30 || peptide.Contains('X'))
            {
                return;
            }
            peptide = peptide.Replace('I', 'L');
            PeptideInfo info = new PeptideInfo { reference = reference, region = region };
            lock (peptideSet)
            {
                if (!peptideSet.ContainsKey(peptide))
                {
                    peptideSet.Add(peptide, new List<PeptideInfo>());
                    peptideSet[peptide].Add(info);
                }
                else
                {
                    List<PeptideInfo> list = peptideSet[peptide];
                    if (list.All((item) => !item.Equals(info)))
                    {
                        peptideSet[peptide].Add(info);
                    }
                }
            }
        }

        static char[] kr = new char[] { 'K', 'R' };
        static int IndexOfKR(string sequence, int startIndex)
        {
            int index = sequence.IndexOfAny(kr, startIndex);
            while (index != -1)
            {
                if (index != sequence.Length - 1 && sequence[index + 1] != 'P')
                {
                    return index;
                }
                else
                {
                    index = sequence.IndexOfAny(kr, index + 1);
                }
            }
            return -1;
        }

        static int IndexOfStopOrKR(string sequence, int startIndex)
        {
            int stopIndex = sequence.IndexOf('*', startIndex);
            int krIndex = IndexOfKR(sequence, startIndex);
            if (stopIndex == -1)
            {
                return krIndex;
            }
            if (krIndex == -1)
            {
                return stopIndex;
            }
            return Math.Min(stopIndex, krIndex);
        }

        static int[] GetSubRegion(int[] source, int start, int len)
        {
            if (source.Length < 2)
            {
                return new int[0];
            }
            start *= 3;
            len *= 3;
            int index = 0;
            int count = source[1] - source[0] + 1;
            while (count <= start)
            {
                start -= count;
                index += 2;
                count = source[index + 1] - source[index] + 1;
            }

            int index2 = index;
            int count2 = count - start;
            while (count2 < len)
            {
                len -= count2;
                index2 += 2;
                count2 = source[index2 + 1] - source[index2] + 1;
            }

            int[] sub = new int[(index2 + 1) - index + 1];
            Array.Copy(source, index, sub, 0, sub.Length);
            sub[0] += start;
            sub[sub.Length - 1] = sub[sub.Length - 2] + len - 1;

            return sub;
        }

        static List<string> tooManyList = new List<string>();
        static int[] shiftList = new int[] { +1, +2, +3, +4, +5, +6, -1, -2, -3, -4, -5, -6, };

        const bool use_predicted_sites = true;

        // return value <染色體編號,ENSG,<peptide序列,座標>>
        static Tuple<string, string, Dictionary<string, HashSet<int[]>>> Work(string genome_number, Strand strand, string reference, Sequence sequence, int[] start, int[] donor, int[] acceptor, int[] ending)
        {
            List<int> donorList = new List<int>(donor);
            foreach (var location in donor)
            {
                if (use_predicted_sites)
                {
                    foreach (var shift in shiftList)
                    {
                        int site = location + shift;
                        if (CheckSpliceSite(sequence, site, 1))
                        {
                            donorList.Add(site);
                        }
                    }
                }
            }
            donor = donorList.Distinct().ToArray();
            List<int> acceptorList = new List<int>(acceptor);
            foreach (var location in acceptor)
            {
                if (use_predicted_sites)
                {
                    foreach (var shift in shiftList)
                    {
                        int site = location + shift;
                        if (CheckSpliceSite(sequence, site, 2))
                        {
                            acceptorList.Add(site);
                        }
                    }
                }
            }
            acceptor = acceptorList.Distinct().ToArray();
            donorList = null;
            acceptorList = null;

            Array.Sort(start);
            Array.Sort(donor);
            Array.Sort(acceptor);

            int count = 0;
            var peptideSet = new Dictionary<string, HashSet<int[]>>();

            sequence.Init();
            var tasks = new Dictionary<int, HashSet<int[]>>();
            foreach (var item in start)
            {
                tasks.Add(item, new HashSet<int[]>(RegionEqualityComparer.Default) { new int[0] });
            }

            //HashSet<int> checkList = new HashSet<int>();

            while (tasks.Count > 0)
            {
                if (count > peptidePerGeneLimit)
                {
                    lock (tooManyList)
                    {
                        tooManyList.Add(reference);
                    }
                    break;
                }
                int startSite = tasks.Keys.Min();
                HashSet<int[]> upstreamSet = tasks[startSite];
                tasks.Remove(startSite);
                //if (upstreamSet.Any((upstream) => upstream.Length == 0) && checkList.Contains(startSite))
                //{
                //    upstreamSet.Remove(new int[0]);
                //}
                //if (upstreamSet.Count == 0)
                //{
                //    continue;
                //}
                foreach (var donorSite in donor.Union(ending).Where((location) => location >= startSite))
                {
                    //foreach (var upstream in upstreamSet)
                    //{
                    //    int frame = CountFrame(upstream);
                    //    char last, next;
                    //    if (frame == 0)
                    //    {
                    //        last = Translation.aminoacids[sequence.SubSequence(acceptorSite, 3)][0];
                    //        next = Translation.aminoacids[sequence.SubSequence(acceptorSite, 3)][0];
                    //    }
                    //    else
                    //    {
                    //        last = Translation.aminoacids[sequence.SubSequence(donorSite - frame + 1, frame) + sequence.SubSequence(acceptorSite, FrameConvert(frame))][0];
                    //        next = Translation.aminoacids[sequence.SubSequence(acceptorSite + FrameConvert(frame), 3)][0];
                    //    }
                    //    if (next != 'P' && (last == 'K' || last == 'R'))
                    //    {
                    //        // lost something
                    //    }
                    //}
                    string aminoAcid = sequence.GetAminoAcid(startSite, donorSite - startSite + 1);
                    int aminoAcidOffset = startSite;
                    bool stopFlag = false;
                    int index = IndexOfStopOrKR(aminoAcid, 0);
                    while (index != -1)
                    {
                        int endSite = index * 3 + aminoAcidOffset + 2;
                        foreach (var upstream in upstreamSet)
                        {
                            int[] region = new int[upstream.Length + 2];
                            Array.Copy(upstream, region, upstream.Length);
                            int frame = CountFrame(upstream);
                            region[region.Length - 2] = startSite - FrameConvert(frame);
                            region[region.Length - 1] = endSite;
                            if (aminoAcid[index] == '*')
                            {
                                region[region.Length - 1] -= 3;
                                if (region[region.Length - 1] < region[region.Length - 2])
                                {
                                    Array.Resize(ref region, region.Length - 2);
                                }
                            }
                            Action<string, int[]> addPeptide = (string target, int[] targetRegion) =>
                            {
                                if (target.Length >= min_peptide_length && target.Length <= max_peptide_length && !target.Contains('X'))
                                {
                                    /*for (int i = 0; i < targetRegion.Length; i += 2)
                                    {
                                        if (targetRegion[i] > targetRegion[i + 1])
                                        {
                                            throw new Exception("座標不正確");
                                        }
                                    }*/
                                    if (strand == Strand.Minus)
                                    {
                                        ReverseLocation(targetRegion, fastaParser.GetSequenceLength(genome_number));
                                    }
                                    if (!peptideSet.ContainsKey(target))
                                    {
                                        peptideSet.Add(target, new HashSet<int[]>(RegionEqualityComparer.Default));
                                    }
                                    peptideSet[target].Add(targetRegion);
                                }
                            };
                            string peptide = GetPeptide(sequence, region).Replace('I', 'L');
                            int krIndex = IndexOfStopOrKR(peptide, 0);
                            if (krIndex == -1)
                            {
                                addPeptide(peptide, region);
                            }
                            else
                            {
                                int startIndex = 0;
                                string fixedPeptide;
                                while (krIndex != -1)
                                {
                                    if (peptide[krIndex] == '*')
                                    {
                                        stopFlag = true;
                                        krIndex -= 1;
                                    }
                                    int len = krIndex - startIndex + 1;
                                    if (len >= min_peptide_length && len <= max_peptide_length)
                                    {
                                        fixedPeptide = peptide.Substring(startIndex, len);
                                        addPeptide(fixedPeptide, GetSubRegion(region, startIndex, len));
                                    }
                                    if (stopFlag)
                                    {
                                        break;
                                    }
                                    startIndex = krIndex + 1;
                                    krIndex = IndexOfStopOrKR(peptide, startIndex);
                                }
                                if (!stopFlag)
                                {
                                    fixedPeptide = peptide.Substring(startIndex);
                                    addPeptide(fixedPeptide, GetSubRegion(region, startIndex, fixedPeptide.Length));
                                }
                            }
                            count++;
                        }
                        if (aminoAcid[index] == '*')
                        {
                            stopFlag = true;
                            break;
                        }
                        upstreamSet.Clear();
                        upstreamSet.Add(new int[0]);
                        startSite = endSite + 1;
                        index = IndexOfStopOrKR(aminoAcid, index + 1);
                    }
                    if (stopFlag)
                    {
                        break;
                    }
                    if (ending.Contains(donorSite))
                    {
                        foreach (var upstream in upstreamSet)
                        {
                            int[] region = new int[upstream.Length + 2];
                            Array.Copy(upstream, region, upstream.Length);
                            int frame = CountFrame(upstream);
                            region[region.Length - 2] = startSite - FrameConvert(frame);
                            region[region.Length - 1] = donorSite;
                            Action<string, int[]> addPeptide = (string target, int[] targetRegion) =>
                            {
                                if (target.Length >= min_peptide_length && target.Length <= max_peptide_length && !target.Contains('X'))
                                {
                                    /*for (int i = 0; i < targetRegion.Length; i += 2)
                                    {
                                        if (targetRegion[i] > targetRegion[i])
                                        {
                                            throw new Exception("座標不正確");
                                        }
                                    }*/
                                    if (strand == Strand.Minus)
                                    {
                                        ReverseLocation(targetRegion, fastaParser.GetSequenceLength(genome_number));
                                    }
                                    if (!peptideSet.ContainsKey(target))
                                    {
                                        peptideSet.Add(target, new HashSet<int[]>(RegionEqualityComparer.Default));
                                    }
                                    peptideSet[target].Add(targetRegion);
                                }
                            };
                            string peptide = GetPeptide(sequence, region).Replace('I', 'L');
                            int krIndex = IndexOfStopOrKR(peptide, 0);
                            if (krIndex == -1)
                            {
                                addPeptide(peptide, region);
                            }
                            else
                            {
                                int startIndex = 0;
                                string fixedPeptide;
                                while (krIndex != -1)
                                {
                                    if (peptide[krIndex] == '*')
                                    {
                                        stopFlag = true;
                                        krIndex -= 1;
                                    }
                                    int len = krIndex - startIndex + 1;
                                    if (len >= min_peptide_length && len <= max_peptide_length)
                                    {
                                        fixedPeptide = peptide.Substring(startIndex, len);
                                        addPeptide(fixedPeptide, GetSubRegion(region, startIndex, len));
                                    }
                                    if (stopFlag)
                                    {
                                        break;
                                    }
                                    startIndex = krIndex + 1;
                                    krIndex = IndexOfStopOrKR(peptide, startIndex);
                                }
                                if (!stopFlag)
                                {
                                    fixedPeptide = peptide.Substring(startIndex);
                                    addPeptide(fixedPeptide, GetSubRegion(region, startIndex, fixedPeptide.Length));
                                }
                            }
                        }
                    }
                    if (!donor.Contains(donorSite))
                    {
                        break;
                    }
                    foreach (var acceptorSite in acceptor.Where((location) => location > donorSite + 1))
                    {
                        foreach (var upstream in upstreamSet)
                        {
                            int[] region = new int[upstream.Length + 2];
                            Array.Copy(upstream, region, upstream.Length);
                            int frame = CountFrame(upstream);
                            region[region.Length - 2] = startSite - FrameConvert(frame);
                            region[region.Length - 1] = donorSite;
                            frame = (frame + region[region.Length - 1] - region[region.Length - 2] + 1) % 3;
                            int fix = acceptorSite + FrameConvert(frame);
                            if (!tasks.ContainsKey(fix))
                            {
                                tasks.Add(fix, new HashSet<int[]>(RegionEqualityComparer.Default) { });
                            }
                            tasks[fix].Add(region);
                        }
                    }
                }
            }
            return new Tuple<string, string, Dictionary<string, HashSet<int[]>>>(genome_number, reference, peptideSet);
        }

        static StreamWriter WriteText(string path)
        {
            const int bufferSize = 8192;
            FileStream fileStream = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.Read, bufferSize, FileOptions.SequentialScan);
            StreamWriter writer = new StreamWriter(fileStream, Encoding.ASCII, bufferSize);
            return writer;
        }

        static void OutputResult(string path, string prefix = "ASG")
        {
            Directory.CreateDirectory(path);
            using (StreamWriter writer = WriteText(Path.Combine(path, "skipList.txt")))
            {
                foreach (string gene in skipList)
                {
                    writer.WriteLine("{0}", gene);
                }
            }
            using (StreamWriter writer = WriteText(Path.Combine(path, "tooManyList.txt")))
            {
                foreach (string gene in tooManyList)
                {
                    writer.WriteLine("{0}", gene);
                }
            }
            using (StreamWriter writer = WriteText(Path.Combine(path, "peptide.fasta")))
            {
                int index = 1;
                foreach (string peptide in peptideSet.Keys)
                {
                    int krIndex = IndexOfKR(peptide, 0);
                    if (krIndex == -1)
                    {
                        writer.WriteLine(">{0}{1:0000000000}", prefix, index);
                        writer.WriteLine(peptide);
                    }
                    else
                    {
                        throw new Exception();
                        //int start = 0;
                        //int fragmentIndex = 1;
                        //while (krIndex != -1)
                        //{
                        //    string fixedPeptide = peptide.Substring(start, fragmentIndex - start + 1);
                        //    if (fixedPeptide.Length >= 5 && fixedPeptide.Length <= 50)
                        //    {
                        //        writer.WriteLine(">{0}{1:0000000000}.{2}", prefix, index, fragmentIndex++);
                        //        writer.WriteLine(fixedPeptide);
                        //    }
                        //    start = fragmentIndex + 1;
                        //    krIndex = IndexOfKR(peptide, krIndex + 1);
                        //}
                    }
                    index++;
                }
            }
            var gzip = new GZipStream(new FileStream(Path.Combine(path, "info.txt.gz"), FileMode.Create, FileAccess.Write), CompressionMode.Compress);
            using (StreamWriter writer = new StreamWriter(gzip))
            {
                int index = 1;
                foreach (var peptide in peptideSet.Values)
                {
                    writer.WriteLine(">{0}{1:0000000000}", prefix, index++);

                    StringBuilder sb = new StringBuilder();

                    foreach (var info in peptide)
                    {
                        var locations = Array.ConvertAll(info.region, (location) => location.ToString());

                        for (int i = 0; i < locations.Length - 1; i++)
                        {
                            sb.Append(locations[i]);
                            if ((i & 0x1) == 0)
                            {
                                sb.Append('-');
                            }
                            else
                            {
                                sb.Append(',');
                            }
                        }
                        sb.Append(locations[locations.Length - 1]);

                        writer.WriteLine("{0}@{1}:{2}", info.reference, info.seq_id, sb.ToString());

                        sb.Clear();
                    }
                }
            }
        }

        static string GetPeptide(Sequence sequence, int[] region)
        {
            string peptide = string.Empty;
            string upstream = string.Empty;
            for (int i = 0; i < region.Length; i += 2)
            {
                if (upstream.Length != 0)
                {
                    string temp = upstream + sequence.GetSequence(region[i], FrameConvert(upstream.Length));
                    if (temp.Contains('N'))
                    {
                        peptide += 'X';
                    }
                    else
                    {
                        string aa = Translation.aminoacids[temp];
                        //if (aa == "*")
                        //{
                        //    break;
                        //}
                        peptide += aa;
                    }
                }
                int start = region[i] + FrameConvert(upstream.Length);
                int length = region[i + 1] - start + 1;
                peptide += sequence.GetAminoAcid(start, length);
                int leave = length % 3;
                upstream = sequence.GetSequence(region[i + 1] - leave + 1, leave);
            }
            return peptide;
        }

        static int FrameConvert(int x)
        {
            switch (x)
            {
                case 0:
                    return 0;
                case 1:
                    return 2;
                case 2:
                    return 1;
                default:
                    throw new Exception();
            }
        }

        static int CountFrame(int[] region)
        {
            int sum = 0;
            for (int i = 0; i < region.Length; i += 2)
            {
                sum += region[i + 1] - region[i] + 1;
            }
            return sum % 3;
        }
    }
}
