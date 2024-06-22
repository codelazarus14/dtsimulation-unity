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
        private bool isRunning;

        private Embed myEmbed;

        // used to force gizmos to only redraw every update 

        void Start()
        {
            MyDT = new DT(config.text);
            StartCoroutine(Run());
        }

        private void OnDrawGizmos()
        {
            int itmp;
            float xPos1, yPos1, xPos2, yPos2;
            // render using gizmos
            Handles.DrawWireDisc(transform.position, transform.forward, scale);

            if (myEmbed == null) return;

            for (int i = 1; i < MyDT.node_number; i++)
            {
                if (myEmbed.Mark[i] != (MyDT.boundary_length + 1))
                {
                    xPos2 = scale * myEmbed.X[i];
                    yPos2 = scale * myEmbed.Y[i];
                    for (int j = MyDT.nstart[i]; j < MyDT.nstart[i + 1]; j++)
                    {
                        itmp = MyDT.ncol[j];
                        xPos1 = scale * myEmbed.X[itmp];
                        yPos1 = scale * myEmbed.Y[itmp];

                        Vector3 start = new Vector2(xPos1, yPos1);
                        Vector3 end = new Vector2(xPos2, yPos2);

                        Handles.DrawAAPolyLine(2f, transform.position + start, transform.position + end);
                    }
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
                    }
                    MyDT.Tidy();
                    MyDT.RelabelNodes();
                    MyDT.Laplacian();

                    myEmbed = new Embed(MyDT);
                    myEmbed.ComputeEmbedding();
                }
                yield return new WaitForSecondsRealtime(timeStep / 1000);
            }
        }
    }
}
