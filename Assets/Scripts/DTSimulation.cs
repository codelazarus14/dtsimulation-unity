using System;
using System.Collections.Generic;
using UnityEngine;

namespace DTSimulation
{
    public class DTSimulation : MonoBehaviour
    {
        // all comments using /**/ are from DT.java

        // TODO: fix style/naming conventions
        // TODO: remove unused

        /* some global constants for run */

        private const float kappa_0 = 0f, BETA = 3f, kappa_d = 0f, kappa_0b = 0f, ALPHA = 0.25f;
        private const int DISKVOL = 472;
        private const int DISKMINVOL = DISKVOL - 30;
        private const int DISKMAXVOL = DISKVOL + 30;
        private const float FRAC = 0.75f;
        private const float TC = 7f;
        private const int MARKEDQ = 282; //FRAC*(DISKVOL+2)/(2.0-FRAC)

        private const int D = 2;
        private const int DPLUS = D + 1;
        private const int DPLUSPLUS = D + 2;
        private const int VOL = DISKVOL + MARKEDQ;
        private const int BIGVOL = 4 * VOL;
        private const int MAXVOL = VOL + 30;
        private const int MINVOL = DPLUSPLUS;
        private const int NUMROW = VOL / 2;
        private const int NONZEROS = NUMROW * VOL; // over 8 billion???

        /* simple global pointers and counters */

        private Node stack_head = null;
        private int simplex_number = 0;
        private int pointer_number = 0;
        private int node_number = 0;

        private int[] boundary;

        private int[] nstart, ncol;
        private double[] nlap;

        /* data structures */

        private Simplex[] simplex_point;
        private int stack_count = 0;

        /* measurements etc */

        private int number_measure = 0, max_point;
        private int[] legal_subsimplex, try_subsimplex, manifold_subsimplex, go_subsimplex;
        private int vol_monitor = 0, num_monitor = 0, b_monitor = 0, growing_vol;

        private bool grow;

        //
        // my fields/deviations from above
        //

        [SerializeField]
        private TextAsset config;


        private void Start() => Main();

        private void Main()
        {
            // abstract graph of points
            // to some coordinate space = enforce smoothness = flatness thru an Action:
            // - maximize dot product of 3d triangle normals? parallel/coplanar
            // - embedded curvature
            // take curr action value, imagine doing it, compute new action resulting, if < old? accept, otherwise accept stochastically... exponential
            // moves towards lower energies with some fluctuations towards higher energy - metropolis algo

            // pick moves at random, link flips or moving vertices
            // generally move towards lower energy w some fluctuations towards higher (depending on temperature)

            int iter, sweep;

            Thermalize();

            //...
        }

        private void Thermalize()
        {
            int therm;

            simplex_point = new Simplex[BIGVOL];
            ReadFile();

            nlap = new double[NONZEROS];
            ncol = new int[NONZEROS];
            nstart = new int[VOL];
            boundary = new int[VOL];

            legal_subsimplex = new int[DPLUS];
            try_subsimplex = new int[DPLUS];
            manifold_subsimplex = new int[DPLUS];
            go_subsimplex = new int[DPLUS];

            /* build lattice */


            //grow = LOGIC.NO;
            //if (grow == LOGIC.YES)
            //{
            //    while (growing_vol < VOL)
            //    {
            //        //System.out.println(growing_vol);
            //        trial_change();
            //        growing_vol += D;
            //        //PrintConfig();
            //    }
            //}

            Tidy();

            //...
        }

        private void ReadFile()
        {
            // TODO: store state in some kind of struct that can be serialized/inspected
            // the temp field, for example, is literally just throwing away some of the Scanner-esque outputs
            // bc those fields are only edited at runtime, not loaded in
            // basically this is all based on a fragile text-based structure that has to be 1-1 with the whitespaced config file
            // when using a structured JSON type thing would save a lot of headaches and look better

            // also TODO: maybe create a textmesh to store these outputs like an in-game console cuz the debug one is not great

            // load config
            string[] configLines = config.text.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
            // using a string stack to imitate use of the Scanner in the java version
            Stack<string> configStack = new Stack<string>();
            for (int i = configLines.Length - 1; i >= 0; i--)
            {
                configStack.Push(configLines[i]);
            }

            int count, s_number, dummy, l_number, c, e1, e2, simp;

            int[] dum = new int[DPLUS];
            int[,] dum2;

            float k0, b;
            float temp;

            // read config values
            try
            {
                s_number = int.Parse(configStack.Pop());
                dum2 = new int[s_number, DPLUS];
                node_number = int.Parse(configStack.Pop());
                count = int.Parse(configStack.Pop());
                // TODO: why do these exist + debug cleanup
                temp = float.Parse(configStack.Pop());
                temp = float.Parse(configStack.Pop());
                temp = float.Parse(configStack.Pop());

                Debug.Log($"{s_number}, {node_number}, {count}");

                simplex_number = 0;
                stack_count = 0;
                pointer_number = 0;

                Debug.Log("Reading in existing configuration");

                Debug.Log($"Configuration has volume: {s_number}");
                Debug.Log($"Node number: {node_number}");
                Debug.Log($"k0b coupling, k2 coupling and curvaturesq coupling: {kappa_0b}, {kappa_d}, {BETA}");

                for (int i = 0; i < count; i++)
                {
                    // TODO: why
                    dummy = int.Parse(configStack.Pop());
                    Push(dummy);
                }

                // set up simplices
                for (int i = 0; i < s_number; i++)
                {
                    for (int j = 0; j < DPLUS; j++)
                    {
                        dum[j] = int.Parse(configStack.Pop());
                        dum2[i, j] = int.Parse(configStack.Pop());
                    }

                    simplex_point[i] = new Simplex(dum, D);
                    simplex_number++;
                    pointer_number++;
                }

                // set neighbors
                for (int i = 0; i < s_number; i++)
                {
                    for (int j = 0; j < DPLUS; j++)
                    {
                        simplex_point[i].neighbors[j] = simplex_point[dum2[i, j]];
                    }
                }

                // print out results
                for (int i = 0; i < s_number; i++)
                {
                    string message = "";
                    for (int j = 0; j < DPLUS; j++)
                    {
                        message += simplex_point[i].vertices[j] + "\t" + simplex_point[i].neighbors[j].label + "\t";
                    }
                    Debug.Log(message);
                }

                growing_vol = simplex_number;

                Debug.Log("Have read data successfully");
            }
            catch (FormatException e)
            {
                Debug.LogError($"Error parsing config file! {e.Message}");
            }
        }

        /* routine pushes deleted vertex label onto stack - DT.java */
        private void Push(int i)
        {
            Node temp = new Node();
            temp.name = i;
            temp.next = stack_head;

            stack_head = temp;
            stack_count++;
            return;
        }

        private void Tidy()
        {
            Simplex[] temp = new Simplex[BIGVOL];
            int add = 0;

            if (pointer_number > max_point) max_point = pointer_number;

            for (int i = 0; i < pointer_number; i++)
            {
                if (simplex_point[i] != null)
                {
                    temp[add] = simplex_point[i];
                    temp[add].label = add;
                    add++;
                }
            }

            for (int i = 0; i < add; i++)
                simplex_point[i] = temp[i];

            pointer_number = add;
            if (pointer_number != simplex_number)
                Debug.Log("oops - pointer number is not equal to simplex_number in tidy()");

            return;
        }
    }
}
