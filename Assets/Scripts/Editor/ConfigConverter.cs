using System;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

namespace DTSimulation.RandomSurface
{
    public class ConfigConverter : AssetPostprocessor
    {
        private const string configsPath = "Assets/Configs";
        private const string outputFolder = "Converted";
        private const string outputPath = configsPath + "/" + outputFolder;

        private static void OnPostprocessAllAssets(string[] importedAssets, string[] deletedAssets, string[] movedAssets, string[] movedFromAssetPaths)
        {
            // after every asset import, check for unconverted config files

            List<(string name, string converted)> converted = new();
            string logMessage = $"-----{nameof(ConfigConverter)}-----\n";

            // check if output folder exists
            if (!AssetDatabase.IsValidFolder(outputPath))
                AssetDatabase.CreateFolder(configsPath, outputFolder);

            // get list of all text assets in the config folder
            string[] configGuids = AssetDatabase.FindAssets("t:textasset", new[] { configsPath });
            foreach (string guid in configGuids)
            {
                string name = GUIDToFilename(guid, false);

                if (!IsConverted(guid))
                {
                    // convert text asset contents into SO
                    string path = AssetDatabase.GUIDToAssetPath(guid);
                    DTConfig config = TextToConfig(AssetDatabase.LoadAssetAtPath<TextAsset>(path));
                    string configPath = $"{outputPath}/{name}.asset";

                    AssetDatabase.CreateAsset(config, configPath);
                    converted.Add((name, $"{name}.asset"));
                }
            }

            foreach (var asset in converted)
                logMessage += $"{asset.name} -> {asset.converted}\n";
            if (converted.Count > 0)
                Debug.Log(logMessage);
        }

        private static DTConfig TextToConfig(TextAsset textAsset)
        {
            DTConfig config = ScriptableObject.CreateInstance<DTConfig>();

            try
            {
                // split string into tokens to be read in one-at-a-time, like DT.ReadFile()
                string[] configTokens = textAsset.text.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                // first two are simplex and node number
                config.simplexCount = int.Parse(configTokens[0]);
                config.nodeCount = int.Parse(configTokens[1]);
                // have to know how many stack lines of the config to skip
                int stackCount = int.Parse(configTokens[2]);

                // skip next 3 (unused coupling params) + number on stack
                int counter = 6 + stackCount;

                config.vertices = new int[config.simplexCount * DT.DPLUS];
                config.neighbors = new int[config.simplexCount * DT.DPLUS];
                // copy simplex data
                for (int i = 0; i < config.simplexCount; i++)
                {
                    for (int j = 0; j < DT.DPLUS; j++)
                    {
                        int configIdx = i * DT.DPLUS + j;
                        config.vertices[configIdx] = int.Parse(configTokens[counter++]);
                        config.neighbors[configIdx] = int.Parse(configTokens[counter++]);
                    }
                }
            }
            catch (FormatException e)
            {
                Debug.LogError($"Error parsing config file! {e.Message}");
            }
            return config;
        }

        private static bool IsConverted(string guid)
        {
            string name = GUIDToFilename(guid, false);
            string[] matches = AssetDatabase.FindAssets($"{name} t:dtconfig", new[] { outputPath });
            // compare filenames in output folder to config's
            foreach (string match in matches)
                if (GUIDToFilename(match, false) == name) return true;
            return false;
        }

        private static string GUIDToFilename(string guid, bool withExtension = true)
        {
            string path = AssetDatabase.GUIDToAssetPath(guid);
            string nameWithExtension = path[(path.LastIndexOf('/') + 1)..];
            return withExtension ? nameWithExtension : nameWithExtension.Remove(nameWithExtension.IndexOf('.'));
        }
    }
}
