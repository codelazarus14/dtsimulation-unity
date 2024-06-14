using UnityEngine;

namespace DTSimulation
{
    public class Simulation : MonoBehaviour
    {
        [SerializeField]
        private bool isRunning;
        [SerializeField]
        private TextAsset config;

        private DT myDT;
        // Start is called before the first frame update
        void Start()
        {
            myDT = new DT(config.text);
            myDT.Thermalize();
        }

        // Update is called once per frame
        void Update()
        {
            if (isRunning)
            {
                for (int i = 0; i < DT.VOL; i++)
                {
                    myDT.TrialChange();
                }
                myDT.Tidy();
                myDT.RelabelNodes();
                myDT.Laplacian();

                //...
            }
        }
    }
}
