using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ASG2
{
    class RegionEqualityComparer : IEqualityComparer<int[]>
    {
        public static RegionEqualityComparer Default { get; } = new RegionEqualityComparer();

        public bool Equals(int[] x, int[] y)
        {
            if (x == null || y == null)
            {
                return x == y;
            }
            return x.SequenceEqual(y);
        }

        public int GetHashCode(int[] obj)
        {
            if (obj == null)
            {
                throw new ArgumentNullException(nameof(obj));
            }
            int hc = obj.Length;
            for (int i = 0, len = hc; i < len; i++)
            {
                unchecked
                {
                    hc = hc * 13 + obj[i];
                }
            }
            return hc;
        }
    }
}
