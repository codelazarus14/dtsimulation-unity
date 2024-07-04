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
        private SerializedProperty wobbleProperty;

        public override VisualElement CreateInspectorGUI()
        {
            VisualElement inspector = inspectorXML.Instantiate();

            betaProperty = serializedObject.FindProperty("beta");
            wobbleProperty = serializedObject.FindProperty("wobbleMagnitude");
            inspector.TrackPropertyValue(betaProperty, OnBetaPropertyChanged);
            inspector.TrackPropertyValue(wobbleProperty, OnWobblePropertyChanged);

            return inspector;
        }

        private void OnBetaPropertyChanged(SerializedProperty property)
        {
            // make beta slider update underlying simulation DT
            DT myDT = (target as Simulation).MyDT;
            if (myDT != null)
                myDT.Beta = property.intValue * 0.1f;
        }

        private void OnWobblePropertyChanged(SerializedProperty property)
        {
            DT myDT = (target as Simulation).MyDT;
            if (myDT != null)
                myDT.WobbleMagnitude = property.floatValue;
        }
    }
}
