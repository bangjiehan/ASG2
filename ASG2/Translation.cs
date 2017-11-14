using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ASG2
{
    class Translation
    {
        public static Dictionary<string, string> aminoacids = new Dictionary<string, string>()
        {
            { "TTT", "F" },
            { "TTC", "F" },
            { "TTA", "L" },
            { "TTG", "L" },
            { "TCT", "S" },
            { "TCC", "S" },
            { "TCA", "S" },
            { "TCG", "S" },
            { "TAT", "Y" },
            { "TAC", "Y" },
            { "TAA", "*" },
            { "TAG", "*" },
            { "TGT", "C" },
            { "TGC", "C" },
            { "TGA", "*" },
            { "TGG", "W" },
            { "CTT", "L" },
            { "CTC", "L" },
            { "CTA", "L" },
            { "CTG", "L" },
            { "CCT", "P" },
            { "CCC", "P" },
            { "CCA", "P" },
            { "CCG", "P" },
            { "CAT", "H" },
            { "CAC", "H" },
            { "CAA", "Q" },
            { "CAG", "Q" },
            { "CGT", "R" },
            { "CGC", "R" },
            { "CGA", "R" },
            { "CGG", "R" },
            { "ATT", "I" },
            { "ATC", "I" },
            { "ATA", "I" },
            { "ATG", "M" },
            { "ACT", "T" },
            { "ACC", "T" },
            { "ACA", "T" },
            { "ACG", "T" },
            { "AAT", "N" },
            { "AAC", "N" },
            { "AAA", "K" },
            { "AAG", "K" },
            { "AGT", "S" },
            { "AGC", "S" },
            { "AGA", "R" },
            { "AGG", "R" },
            { "GTT", "V" },
            { "GTC", "V" },
            { "GTA", "V" },
            { "GTG", "V" },
            { "GCT", "A" },
            { "GCC", "A" },
            { "GCA", "A" },
            { "GCG", "A" },
            { "GAT", "D" },
            { "GAC", "D" },
            { "GAA", "E" },
            { "GAG", "E" },
            { "GGT", "G" },
            { "GGC", "G" },
            { "GGA", "G" },
            { "GGG", "G" },
        };

        public static string Translate(string sequence, int start)
        {
            StringBuilder sb = new StringBuilder();
            int index = start == -1 ? sequence.IndexOf("ATG") : start;
            string aminoacid;
            do
            {
                string codon = sequence.Substring(index, 3);
                if (codon.Contains('N'))
                {
                    aminoacid = "X";
                }
                else
                {
                    aminoacid = aminoacids[codon];
                }
                sb.Append(aminoacid);
                index += 3;
            } while (index + 2 <= sequence.Length - 1);
            return sb.ToString();
        }

        public static char TranslateChar(string sequence)
        {
            if (sequence.Contains('N'))
            {
                return 'X';
            }
            else
            {
                return aminoacids[sequence][0];
            }
        }

        public static string Translate(string sequence)
        {
            return Translate(sequence, -1);
        }

        public static IEnumerable<string> KR(string aminoAcidSequence)
        {
            int index = 0;
            int index2 = 0;

            while (true)
            {
                index2 = indexOfKR(aminoAcidSequence, index);
                if (index2 == -1)
                {
                    if (index != aminoAcidSequence.Length)
                    {
                        yield return aminoAcidSequence.Substring(index);
                    }
                    break;
                }
                yield return aminoAcidSequence.Substring(index, index2 - index + 1);
                index = index2 + 1;
            }
        }

        public static int indexOfKR(string aminoAcidSequence, int startIndex)
        {
            for (int i = startIndex; i < aminoAcidSequence.Length; i++)
            {
                if ((aminoAcidSequence[i] == 'K' || aminoAcidSequence[i] == 'R') && (i == aminoAcidSequence.Length - 1 || aminoAcidSequence[i + 1] != 'P'))
                {
                    return i;
                }
            }
            return -1;
        }
    }
}
