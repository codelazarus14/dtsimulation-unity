using System;
using System.Collections;
using UnityEditor;
using UnityEngine;

namespace DTSimulation.RandomSurface
{
    public class Simulation : MonoBehaviour
    {
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
        [Range(0.01f, 1f)]
        private float nodeSize = 0.1f;
        [SerializeField]
        private bool trialChange;
        [SerializeField]
        private bool drawAveragePos;
        [SerializeField]
        private bool drawNormals;
        [SerializeField]
        private bool isRunning;

        void Start()
        {
            MyDT = new DT(config.text);
            StartCoroutine(Run());
        }

        private void OnDrawGizmos()
        {
            if (MyDT == null) return;

            Vector3 averagePos = Vector3.zero;

            for (int i = 0; i < MyDT.node_number; i++)
            {
                Vector3 node2Pos = scale * MyDT.NodePositions[i];
                Gizmos.color = Color.red;
                Gizmos.DrawSphere(node2Pos, HandleUtility.GetHandleSize(node2Pos) * nodeSize);

                averagePos += node2Pos;

                (int[] neighbors, int nCount) = MyDT.NearestNeighbors(i);

                for (int j = 0; j < nCount; j++)
                {
                    int neighbor = neighbors[j];
                    Vector3 node1Pos = scale * MyDT.NodePositions[neighbor];

                    Handles.DrawAAPolyLine(2f, transform.position + node1Pos, transform.position + node2Pos);
                }
            }

            // draw normals
            if (drawNormals)
            {
                for (int i = 0; i < MyDT.pointer_number; i++)
                {
                    Simplex p = MyDT.simplex_point[i];
                    if (p == null) continue;

                    // compute normal from verts
                    Vector3 center = Vector3.zero;
                    foreach (int v in p.vertices)
                        center += MyDT.NodePositions[v];
                    center /= p.vertices.Length;

                    // TODO: fix (this doesn't always keep the same orientation across triangles)
                    Vector3 n = Vector3.Cross(
                        MyDT.NodePositions[p.vertices[1]] - MyDT.NodePositions[p.vertices[0]],
                        MyDT.NodePositions[p.vertices[2]] - MyDT.NodePositions[p.vertices[0]]).normalized;
                    Vector3 normalPos = (transform.position + center) * scale;
                    Vector3 normalTip = normalPos + n * HandleUtility.GetHandleSize(normalPos) * 2f;

                    Handles.color = Color.blue;
                    Handles.DrawAAPolyLine(5f, normalPos, normalTip);
                    Gizmos.color = Color.white;
                    Gizmos.DrawSphere(normalTip, HandleUtility.GetHandleSize(normalPos) * 0.05f);
                    Handles.color = Color.white;
                }
            }

            averagePos /= MyDT.node_number;

            if (drawAveragePos)
            {
                Gizmos.color = Color.green;
                Gizmos.DrawSphere(averagePos, HandleUtility.GetHandleSize(averagePos));
                Gizmos.color = Color.white;
            }
        }

        private IEnumerator Run()
        {
            MyDT.Thermalize();
            MyDT.RelabelNodes();
            while (true)
            {
                if (isRunning)
                {
                    // main work of simulation
                    for (int i = 0; i < DT.VOL; i++)
                    {
                        if (trialChange)
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
