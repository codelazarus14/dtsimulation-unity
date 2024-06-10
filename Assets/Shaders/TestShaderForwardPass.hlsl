
// used to get basic functions like GetVertexPositionInputs
#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"

// mesh data passed into the shader
struct Attributes
{
    float3 position : POSITION;
};

struct Interpolators
{
    float4 positionCS : SV_POSITION;
};

CBUFFER_START(UnityPerMaterial)
    half4 _BaseColor;
CBUFFER_END

// vertex shader
Interpolators Vertex(Attributes input)
{
    Interpolators output;
    
    VertexPositionInputs positions = GetVertexPositionInputs(input.position);
    
    float4 clipSpacePosn = positions.positionCS;
    output.positionCS = clipSpacePosn;
    
    return output;
}

// fragment shader
float4 Fragment(Interpolators input) : SV_TARGET
{
    return _BaseColor;
}