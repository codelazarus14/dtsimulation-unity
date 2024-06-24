using System;
using System.Collections;
using UnityEditor;
using UnityEngine;

namespace DTSimulation.RandomSurface
{
    public class Simulation : MonoBehaviour
    {
        // TODO: make new editor script to have this one work as well
        public DT MyDT { get; private set; }

        [SerializeField]
        private TextAsset config;
        // speed in msec
        [SerializeField]
        [Range(50f, 500f)]
        private float timeStep = 100f;
        [SerializeField]
        [Range(0, 100)]
        private int beta = 30;
        [SerializeField]
        [Range(1f, 10f)]
        private float scale = 5f;
        [SerializeField]
        private bool isRunning;


        void Start()
        {
            MyDT = new DT(config.text);
            StartCoroutine(Run());
        }

        private void OnDrawGizmos()
        {
            Gizmos.DrawWireSphere(transform.position, scale);
            if (MyDT == null) return;

            for (int i = 1; i < MyDT.node_number; i++)
            {
                Vector3 node2Pos = scale * MyDT.NodePositions[i];
                Gizmos.color = Color.red;
                // TODO: node radius parameter
                Gizmos.DrawSphere(node2Pos, 0.1f);
                Gizmos.color = Color.white;

                (int[] neighbors, int nCount) = MyDT.GetNodeNN(i);

                for (int j = 0; j < nCount; j++)
                {
                    int neighbor = neighbors[j];
                    Vector3 node1Pos = scale * MyDT.NodePositions[neighbor];

                    Handles.DrawAAPolyLine(2f, transform.position + node1Pos, transform.position + node2Pos);
                }
            }
        }

        private IEnumerator Run()
        {
            MyDT.Thermalize();
            while (true)
            {
                if (isRunning)
                {
                    // main work of simulation
                    for (int i = 0; i < DT.VOL; i++)
                    {
                        MyDT.TrialChange();
                        // wobble vertex positions
                        MyDT.WobbleVertex();
                    }
                    MyDT.Tidy();
                }
                yield return new WaitForSecondsRealtime(timeStep / 1000);
            }
        }
    }
}
