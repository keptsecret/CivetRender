#version 420 core
out vec4 FragColor;

uniform sampler2D PositionMap;
uniform sampler2D AlbedoMap;
uniform sampler2D AORoughMetallicMap;
uniform sampler2D NormalMap;
uniform vec2 screenSize;
uniform vec3 viewPos;

uniform ivec3 probeGridDims;
uniform vec3 gridCornerCoord;
uniform vec3 gridCellSize;
uniform int numIrradianceSamples;
uniform int numDistanceSamples;
uniform float irradianceLobeSize;
uniform float distanceLobeSize;
uniform sampler2DArray radianceOctMap;
uniform sampler2DArray distanceOctMap;

const float Pi = 3.141592654f;

int gridCoordToFlatIndex(vec3 probePos) {
    return int(probePos.x + probeGridDims.x * probePos.y + probeGridDims.x * probeGridDims.y * probePos.z);
}

ivec3 baseGridCoord(vec3 pos) {
    return clamp(ivec3((pos - gridCornerCoord) / gridCellSize),
            ivec3(0),
            ivec3(probeGridDims) - ivec3(1));
}

int baseProbeIndex(vec3 pos) {
    return gridCoordToFlatIndex(baseGridCoord(pos));
}

ivec3 flatIndexToGridCoord(int index) {
    return ivec3(mod(index, probeGridDims.x),
    mod(index / probeGridDims.x, probeGridDims.y),
    index / (probeGridDims.x * probeGridDims.y));
}

vec3 gridCoordToPosition(ivec3 coord) {
    return gridCellSize * vec3(coord) + gridCornerCoord;
}

vec3 getProbePosition(int index) {
    return gridCoordToPosition(flatIndexToGridCoord(index));
}

float signNotZero(float f) {
    return(f >= 0.0) ? 1.0 : -1.0;
}
vec2 signNotZero(vec2 v) {
    return vec2(signNotZero(v.x), signNotZero(v.y));
}

vec2 octEncode(vec3 v) {
    float l1norm = abs(v.x) + abs(v.y) + abs(v.z);
    vec2 result = v.xy * (1.0 / l1norm);
    if (v.z < 0.0) {
        result = (1.0 - abs(result.yx)) * signNotZero(result.xy);
    }
    return result;
}

vec3 octDecode(vec2 o) {
    vec3 v = vec3(o.x, o.y, 1.0 - abs(o.x) - abs(o.y));
    if (v.z < 0.0) {
        v.xy = (1.0 - abs(v.yx)) * signNotZero(v.xy);
    }
    return normalize(v);
}

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

vec3 fresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness) {
    return F0 + (max(vec3(1.0 - roughness), F0) - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

void main() {
    vec2 texCoords = getTexCoords();
    vec3 albedo = texture(AlbedoMap, texCoords).xyz;
    float ao = texture(AORoughMetallicMap, texCoords).r;
    float roughness = texture(AORoughMetallicMap, texCoords).g;
    float metallic = texture(AORoughMetallicMap, texCoords).b;

    vec3 worldPos = texture(PositionMap, texCoords).xyz;
    vec3 normal = normalize(texture(NormalMap, texCoords).xyz);
    vec3 viewDir = normalize(viewPos - worldPos);

    vec3 F0 = vec3(0.04);
    F0 = mix(F0, albedo, metallic);
    vec3 kS = fresnelSchlickRoughness(max(dot(normal, viewDir), 0.0), F0, roughness);
    vec3 kD = 1.0 - kS;
    kD *= 1.0 - metallic;
    //vec3 specularAlbedo = mix(vec3(0.003), albedo, vec3(metallic));

    ivec3 baseCoord = baseGridCoord(worldPos);
    vec3 basePos = gridCoordToPosition(baseCoord);
    vec3 alpha = clamp((worldPos - basePos) / gridCellSize, vec3(0.0), vec3(1.0));  // interpolation weights

    // selects 8 probes surrounding fragment
    vec3 irradiance = vec3(0.0);
    float sum_weight = 0.0;
    for (int idx = 0; idx < 8; idx++) {
        ivec3 offset = ivec3(idx, idx >> 1, idx >> 2) & ivec3(1);
        ivec3 probeGridCoord = ivec3(clamp(baseCoord + offset, vec3(0.0), vec3(probeGridDims.x, probeGridDims.y, probeGridDims.z) - vec3(1.0)));
        int probeIndex = gridCoordToFlatIndex(probeGridCoord);

        vec3 trilinear = mix(1 - alpha, alpha, offset);
        float weight = trilinear.x * trilinear.y * trilinear.z;

        vec3 probePos = gridCoordToPosition(probeGridCoord);
        vec3 probeDir = normalize(worldPos - probePos);

        weight *= max(0.05, dot(-probeDir, normal));

        vec2 temp = texture(distanceOctMap, vec3(octEncode(probeDir) * 0.5 + 0.5, probeIndex)).xy;
        float mean = temp.x;
        float variance = abs(temp.y - (mean * mean));

        float probeDist = length(worldPos - probePos);
        float t_sub_mean = probeDist - mean;
        float chebychev = variance / (variance + (t_sub_mean * t_sub_mean));
        weight *= ((probeDist <= mean)) ? 1.0 : max(0.0, chebychev);

        weight = max(0.002, weight);
        sum_weight += weight;

        vec3 indirectDiffuse = texture(radianceOctMap, vec3(octEncode(normal) * 0.5 + 0.5, probeIndex)).xyz;

        irradiance += indirectDiffuse * albedo * kD * weight;
    }

    FragColor.rgb = 2.0 * Pi * irradiance * ao / sum_weight;
}

