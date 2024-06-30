using UnityEditor;
using UnityEditor.UIElements;
using UnityEngine.UIElements;

namespace DTSimulation.RandomSurface
{
    [CustomEditor(typeof(Simulation))]
    public class SimulationRS_Inspector : Editor
    {
        public VisualTreeAsset inspectorXML;

        private SerializedProperty betaProperty;

        public override VisualElement CreateInspectorGUI()
        {
            VisualElement inspector = inspectorXML.Instantiate();

            betaProperty = serializedObject.FindProperty("beta");
            inspector.TrackPropertyValue(betaProperty, OnBetaPropertyChanged);

            return inspector;
        }

        private void OnBetaPropertyChanged(SerializedProperty property)
        {
            // make beta slider update underlying simulation DT
            DT myDT = (target as Simulation).MyDT;
            if (myDT != null)
                myDT.Beta = property.intValue * 0.1f;
        }
    }
}
