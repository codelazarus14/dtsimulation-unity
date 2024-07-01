using UnityEngine;

namespace DTSimulation.RandomSurface
{
    public class DTConfig : ScriptableObject
    {
        public int simplexCount;
        public int nodeCount;
        public int[] vertices;
        public int[] neighbors;
    }
}
