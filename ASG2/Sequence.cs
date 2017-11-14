using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ASG2
{
    class Sequence
    {
        string sequence;
        int offset;
        string[] aminoAcid;

        public int Start
        {
            get
            {
                return offset + 1;
            }
        }

        public int End
        {
            get
            {
                return sequence.Length + offset;
            }
        }

        public int Length
        {
            get
            {
                return sequence.Length;
            }
        }

        public Sequence(string sequence, int offset)
        {
            this.sequence = sequence;
            this.offset = offset;
        }

        public bool Init()
        {
            if (aminoAcid == null)
            {
                aminoAcid = new string[3];
                for (int i = 0; i < 3; i++)
                {
                    aminoAcid[i] = Translation.Translate(sequence, i);
                }
                return true;
            }
            return false;
        }

        public string GetSequence(int start, int length)
        {
            start -= offset;
            return sequence.Substring(start - 1, length);
        }

        public string GetAminoAcid(int start, int length)
        {
            start -= offset;
            int frame;
            int index = Math.DivRem(start - 1, 3, out frame);
            return aminoAcid[frame].Substring(index, length / 3);
        }
    }
}
