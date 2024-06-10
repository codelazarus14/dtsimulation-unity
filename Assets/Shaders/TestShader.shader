Shader "Unlit/TestShader"
{
	Properties 
	{
		_BaseColor("Base Color", Color) = (1, 1, 1, 1)
	}
	SubShader 
	{
		// defines this as the shader for the universal pipeline
		Tags{"RenderPipeline" = "UniversalPipeline"}

		Pass
		{
			// this is the main lighting pass, of this subshader
			Name "ForwardPass"
			Tags{"LightMode" = "UniversalForward"}

			HLSLPROGRAM

			#pragma vertex Vertex
			#pragma fragment Fragment

			#include "TestShaderForwardPass.hlsl"

			ENDHLSL
		}
	}
}
