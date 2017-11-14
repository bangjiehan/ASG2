using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ASG2
{
    public enum Strand : byte
    {
        None,
        Plus,
        Minus,
        Unknown,
    }

    class GffParser
    {
        public class Node
        {
            public string seqid;
            public string source;
            public string feature;
            public int start;
            public int end;
            public float? score;
            public Strand strand;
            public byte? phase;

            // GFF column 9 attributes table, null
            public Dictionary<string, string> attributes;

            // 首字大寫皆為保留屬性, 首字小寫可自由使用
            // ID只僅允許一個值, Parent Alias Note Dbxref Ontology_term允許多個值

            // Just for convenient access.
            public string ID
            {
                get
                {
                    return attributes.ContainsKey(GffParser.ID) ? attributes[GffParser.ID] : null;
                }
            }

            // 多個值以逗號分隔
            public string Parents
            {
                get
                {
                    return attributes.ContainsKey(PARENT) ? attributes[PARENT] : null;
                }
            }

            // GFF資料行
            public string GffLineString
            {
                get
                {
                    StringBuilder buffer = new StringBuilder();
                    char strandChar;
                    switch (strand)
                    {
                        case Strand.None:
                            strandChar = '.';
                            break;
                        case Strand.Plus:
                            strandChar = '+';
                            break;
                        case Strand.Minus:
                            strandChar = '-';
                            break;
                        case Strand.Unknown:
                            strandChar = '?';
                            break;
                        default:
                            strandChar = '.';
                            break;
                    }

                    buffer.AppendFormat("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t", seqid ?? ".", source ?? ".", feature ?? ".", start, end, score.HasValue ? score.ToString() : ".", strandChar, phase.HasValue ? phase.ToString() : ".");

                    foreach (var attribute in attributes)
                    {
                        buffer.Append(attribute.Key + "=" + attribute.Value + ";");
                    }

                    return buffer.ToString();
                }
            }

            // Sub Nodes
            public List<Node> subNodes;
            // Parent Nodes
            public List<Node> parentNodes;

            // Init
            public Node()
            {
                attributes = new Dictionary<string, string>();
                subNodes = new List<Node>();
                parentNodes = new List<Node>();
            }

            public Node GetSubNode(string id)
            {
                return subNodes.Find((node) => node.ID == id);
            }
        }

        // For saving memory.
        // seqid, source, type table

        private const string GZ = ".gz";
        private const string ID = "ID";
        private const string PARENT = "Parent";

        // Exclude list
        List<string> excludeFeatures;
        List<string> attributeFilter;

        // Strand lookup table
        readonly static ReadOnlyDictionary<string, Strand> strandLookupTable = new ReadOnlyDictionary<string, Strand>(new Dictionary<string, Strand>()
        {
            { ".", Strand.None },
            { "+", Strand.Plus },
            { "-", Strand.Minus },
            { "?", Strand.Unknown },
        });

        public GffParser()
        {
            excludeFeatures = new List<string> { "exon", "CDS", "UTR", "five_prime_UTR", "three_prime_UTR", };
            attributeFilter = new List<string> { "ID", "Parent", };
        }

        public GffParser(string[] attributeFilter) : this()
        {
            if (attributeFilter != null)
            {
                this.attributeFilter.AddRange(attributeFilter);
                this.attributeFilter = this.attributeFilter.Distinct().ToList();
            }
        }

        public GffParser(string[] attributeFilter, string[] excludeArray) : this(attributeFilter)
        {
            if (excludeArray != null)
            {
                excludeFeatures = excludeArray.Distinct().ToList();
            }
            else
            {
                excludeFeatures = null;
            }
        }

        private Node CreateNode(string line)
        {
            Node node = new Node();

            const int ColumnLength = 9;
            string[] columns = line.Split(new char[] { '\t' }, ColumnLength + 1);
            if (columns.Length != ColumnLength)
            {
                throw new Exception("GFF欄數錯誤");
            }

            const string DOT = ".";
            if (columns[0] != DOT)
            {
                node.seqid = columns[0];
            }
            if (columns[1] != DOT)
            {
                node.source = columns[1];
            }
            if (columns[2] != DOT)
            {
                node.feature = columns[2];
            }
            if (columns[3] != DOT)
            {
                node.start = int.Parse(columns[3]);
            }
            if (columns[4] != DOT)
            {
                node.end = int.Parse(columns[4]);
            }
            if (columns[5] != DOT)
            {
                node.score = float.Parse(columns[5]);
            }
            if (columns[6] != DOT)
            {
                node.strand = strandLookupTable[columns[6]];
            }
            if (columns[7] != DOT)
            {
                node.phase = byte.Parse(columns[7]);
            }

            // Attribute table
            if (columns[8] != DOT)
            {
                foreach (string attribute in columns[8].Split(new char[] { ';' }, StringSplitOptions.RemoveEmptyEntries))
                {
                    string[] pair = attribute.Split(new char[] { '=' }, 2 + 1);
                    if (pair.Length != 2)
                    {
                        // 與GFF2(GTF)相容
                        //pair = attribute.Split(new char[] { ' ' }, 2 + 1, StringSplitOptions.RemoveEmptyEntries);
                        //if (pair.Length != 2)
                        //{
                        throw new Exception("屬性欄數異常");
                        //}
                    }
                    // 過濾屬性
                    if (attributeFilter != null && !attributeFilter.Contains(pair[0]))
                    {
                        continue;
                    }

                    node.attributes.Add(pair[0], pair[1]);
                }
            }

            // 檢查ID屬性
            if (node.attributes.ContainsKey(ID) && node.attributes[ID].IndexOf(',') != -1)
            {
                throw new Exception("ID屬性含有多重值");
            }

            return node;
        }

        public string[] ExcludeFeature
        {
            get
            {
                return excludeFeatures.ToArray();
            }
        }

        public string[] AttributeFilter
        {
            get
            {
                return attributeFilter.ToArray();
            }
        }

        public List<Node> ReadAllNodes(string gffFilePath)
        {
            Stream stream = File.OpenRead(gffFilePath);

            if (gffFilePath.EndsWith(GZ))
            {
                stream = new GZipStream(stream, CompressionMode.Decompress);
            }

            return ReadAllNodes(stream);
        }

        public List<Node> ReadAllNodes(Stream gffStream)
        {
            List<Task> tasks = new List<Task>();

            using (var reader = new StreamReader(gffStream, Encoding.ASCII, false, 4096))
            {
                // GFF3版本
                string[] version = reader.ReadLine().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                if (version[0] != "##gff-version" || !version[1].StartsWith("3"))
                {
                    // GTF沒有此標記
                    throw new Exception("GFF3版本錯誤");
                }

                while (!reader.EndOfStream)
                {
                    List<string> lineList = new List<string>();
                    string line;
                    do
                    {
                        line = reader.ReadLine();
                        if (line.StartsWith("#"))
                        {
                            if (line == "###") // ### 代表feature的屬性都不會再參照之後的feature
                            {
                                //if (lineList.Count >= 1000) // 至少有1000筆才建立新工作
                                {
                                    break;
                                }
                            }
                            else
                            {
                                continue;
                            }
                        }
                        else
                        {
                            lineList.Add(line);
                        }
                    } while (!reader.EndOfStream);

                    // 多執行緒處理
                    tasks.Add(Task.Run(() =>
                    {
                        return ProcessNodes(lineList);
                    }));
                }
            }
            Task.WaitAll(tasks.ToArray());

            List<Node> nodes = new List<Node>();

            foreach (Task<IEnumerable<Node>> task in tasks)
            {
                nodes.AddRange(task.Result);
            }

            return nodes;
        }

        private IEnumerable<Node> ProcessNodes(IEnumerable<string> lineList)
        {
            // Nodes contained ID
            Dictionary<string, Node> nodeTable = new Dictionary<string, Node>();
            // All nodes
            List<Node> allNodes = new List<Node>();

            // Root nodes
            List<Node> rootNodes = new List<Node>();

            foreach (string line in lineList)
            {
                Node node = CreateNode(line);
                if (node.ID != null && (excludeFeatures == null || !excludeFeatures.Contains(node.feature))) // 排除包含多重節點的CDS,假設CDS沒有任何子節點
                {
                    nodeTable.Add(node.ID, node);
                }
                allNodes.Add(node);
            }

            // set node relationship
            foreach (Node node in allNodes)
            {
                if (node.Parents == null)
                {
                    if (node.ID == null)
                    {
                        throw new Exception("Root節點沒有ID");
                    }
                    rootNodes.Add(node);
                }
                else
                {
                    if (node.Parents.IndexOf(',') == -1) // not contain ','
                    {
                        LinkNodes(nodeTable[node.Parents], node);
                    }
                    else
                    {
                        throw new Exception("不允許多個父節點");
                    }
                }
            }
            return rootNodes;
            //return allNodes;
        }

        // Link nodes.
        private static void LinkNodes(Node parentNode, Node subNode)
        {
            parentNode.subNodes.Add(subNode);
            subNode.parentNodes.Add(parentNode);
        }
    }
}
