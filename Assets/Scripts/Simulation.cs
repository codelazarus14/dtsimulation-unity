using System;
using System.Collections;
using UnityEditor;
using UnityEngine;

namespace DTSimulation
{
    public class Simulation : MonoBehaviour
    {
        [SerializeField]
        private bool isRunning;
        [SerializeField]
        private TextAsset config;
        // speed in msec
        [SerializeField]
        [Range(50f, 500f)]
        private float speed = 200f;
        [SerializeField]
        [Range(1f, 300f)]
        private float scale;
        // TODO: add DT.BETA slider

        private DT myDT;
        private Embed myEmbed;

        // used to force gizmos to only redraw every update 

        void Start()
        {
            myDT = new DT(config.text);
            StartCoroutine(Run());
        }

        private void OnDrawGizmos()
        {
            int itmp;
            float xPos1, yPos1, xPos2, yPos2;
            // render using gizmos
            Handles.DrawWireDisc(transform.position, transform.forward, scale);

            if (myEmbed == null) return;

            for (int i = 1; i < myDT.node_number; i++)
            {
                if (myEmbed.Mark[i] != (myDT.boundary_length + 1))
                {
                    xPos2 = scale * myEmbed.X[i];
                    yPos2 = scale * myEmbed.Y[i];
                    for (int j = myDT.nstart[i]; j < myDT.nstart[i + 1]; j++)
                    {
                        itmp = myDT.ncol[j];
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
            myDT.Thermalize();
            while (true)
            {
                if (isRunning)
                {
                    // main work of simulation
                    for (int i = 0; i < DT.VOL; i++)
                    {
                        myDT.TrialChange();
                    }
                    myDT.Tidy();
                    myDT.RelabelNodes();
                    myDT.Laplacian();

                    myEmbed = new Embed(myDT);
                    myEmbed.ComputeEmbedding();
                }
                yield return new WaitForSecondsRealtime(speed / 1000);
            }
        }
    }
}
