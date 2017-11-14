/*
 * FastaParser 2.1 for ASG2
 * 2016/6/13
 * Bang-Jie, Han
 */

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ASG2
{
    class FastaParser
    {
        struct Index
        {
            public int Length;
            public long Offset;
            public int LineBases;
            public int LineWidth;
        }

        StreamReader fastaReader;
        Dictionary<string, Index> fastaIndex;

        public FastaParser(string path) : this(path, false)
        {
        }

        public FastaParser(string path, bool isAsync)
        {
            fastaReader = OpenText(path, (isAsync ? FileOptions.Asynchronous | FileOptions.RandomAccess : FileOptions.RandomAccess));
            fastaIndex = ParseIndex(Path.GetFullPath(path) + ".fai");
        }

        static Dictionary<string, Index> ParseIndex(string path)
        {
            var fastaIndex = new Dictionary<string, Index>();
            using (var reader = OpenText(path, FileOptions.SequentialScan))
            {
                while (!reader.EndOfStream)
                {
                    string[] line = reader.ReadLine().Split('\t');
                    Index index;
                    index.Length = int.Parse(line[1]);
                    index.Offset = long.Parse(line[2]);
                    index.LineBases = int.Parse(line[3]);
                    index.LineWidth = int.Parse(line[4]);
                    fastaIndex.Add(line[0], index);
                }
            }
            return fastaIndex;
        }

        static StreamReader OpenText(string path, FileOptions options)
        {
            const int bufferSize = 4096;
            var file = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, bufferSize, options);
            var reader = new StreamReader(file, Encoding.ASCII, false, bufferSize, false);
            return reader;
        }

        public IEnumerable<string> GetSequenceNames()
        {
            return fastaIndex.Keys;
        }

        public int GetSequenceLength(string sequenceName)
        {
            return fastaIndex[sequenceName].Length;
        }

        public string Read(string sequenceName, int start, int count)
        {
            return Read(sequenceName, start, count, false);
        }

        public string Read(string sequenceName, int start, int count, bool reverse)
        {
            char[] buffer = ReadInternal(sequenceName, start - 1, count);
            if (reverse)
            {
                Reverse(buffer);
            }
            return new string(buffer);
        }

        public string ReadBlock(string sequenceName, int start, int end)
        {
            return ReadBlock(sequenceName, start, end, false);
        }

        public string ReadBlock(string sequenceName, int start, int end, bool reverse)
        {
            int offset = start - 1;
            char[] buffer = ReadInternal(sequenceName, offset, end - offset);
            if (reverse)
            {
                Reverse(buffer);
            }
            return new string(buffer);
        }

        char[] ReadInternal(string sequenceName, int offset, int count)
        {
            Index index = fastaIndex[sequenceName];
            int lineBases = index.LineBases;
            int offsetInLine;
            int lineOffset = Math.DivRem(offset, lineBases, out offsetInLine);
            fastaReader.BaseStream.Position = index.Offset + lineOffset * index.LineWidth + offsetInLine;
            fastaReader.DiscardBufferedData();
            char[] buffer = new char[count];
            int firstLineBases = lineBases - offsetInLine;
            if (count <= firstLineBases)
            {
                fastaReader.Read(buffer, 0, count);
            }
            else
            {
                int bufferindex = fastaReader.Read(buffer, 0, firstLineBases);
                int newlineWidth = index.LineWidth - lineBases;
                char[] nullBuffer = new char[newlineWidth];
                fastaReader.Read(nullBuffer, 0, newlineWidth);
                while (bufferindex + lineBases < count)
                {
                    bufferindex += fastaReader.Read(buffer, bufferindex, lineBases);
                    fastaReader.Read(nullBuffer, 0, newlineWidth);
                }
                fastaReader.Read(buffer, bufferindex, count - bufferindex);
            }
            return buffer;
        }

        static void Reverse(char[] buffer)
        {
            Array.Reverse(buffer);
            for (int i = 0, len = buffer.Length; i < len; i++)
            {
                switch (buffer[i])
                {
                    case 'A':
                        buffer[i] = 'T';
                        break;
                    case 'T':
                        buffer[i] = 'A';
                        break;
                    case 'C':
                        buffer[i] = 'G';
                        break;
                    case 'G':
                        buffer[i] = 'C';
                        break;
                }
            }
        }
    }
}
