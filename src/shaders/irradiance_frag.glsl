#version 420 core
layout (location = 0) out vec4 irradianceOut;

const int NUM_SPHERE_SAMPLES = 128;
layout (std140, binding = 2) uniform SphereDirectionSamples {
    vec3 sphereSamples[NUM_SPHERE_SAMPLES];
};

in vec2 TexCoord;

uniform samplerCube cubemap;
uniform int probeIndex;
uniform int numSamples;
uniform float lobeSize;

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

vec3 directionInSphere(int i, int n) {
    int k = NUM_SPHERE_SAMPLES / n;

    int index = i * k;
    vec3 dir = sphereSamples[index % NUM_SPHERE_SAMPLES].xyz;
    return dir;
}

void main() {
    vec3 baseDir = octDecode(TexCoord * 2.0 - 1.0);
    vec3 irradiance = vec3(0.0);

    for (int i = 0; i < numSamples; i++) {
        vec3 offset = directionInSphere(i, numSamples);
        vec3 sampleDir = normalize(baseDir + lobeSize * offset);
        vec3 sampleRadiance = texture(cubemap, sampleDir).rgb;

        irradiance += sampleRadiance;
    }

    irradiance /= float(numSamples);
    irradianceOut = vec4(irradiance, 1.0);
}