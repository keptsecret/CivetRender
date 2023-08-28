#version 420 core
out vec4 FragColor;

uniform sampler2D DepthMap;
uniform sampler2D NormalMap;
uniform sampler2D AORoughMetallicMap;
uniform sampler2D RawFinalImage;
uniform samplerCube skyboxSampler;

uniform vec2 screenSize;
uniform mat4 projToPixel;
uniform mat4 sceneView;
uniform mat4 sceneInvView;
uniform mat4 sceneProjection;
uniform mat4 sceneInvProjection;
uniform float nearPlane;
uniform float farPlane;

// SSR parameters
const float maxSteps = 50;
const float binarySearchIterations = 5;
const float jitterAmount = 1.0;
const float maxDistance = 200.0;
const float stride = 8.0;
const float zThickness = 1.5;
const float strideZCutoff = 100.0;
const float screenEdgeFadeStart = 0.75;
const float eyeFadeStart = 0.5;
const float eyeFadeEnd = 1.0;

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

float distanceSquared(vec2 a, vec2 b) {
    a -= b;
    return dot(a, a);
}

float linear01Depth(float z) {
    float tmp = farPlane / nearPlane;
    return 1.0 / ((1.0 - tmp) * z + tmp);
}

// adapted from McGuire's DDA solution: https://casual-effects.blogspot.com/2014/08/screen-space-ray-tracing.html
bool raymarch(vec3 csOrig, vec3 csDir, float jitter, out vec2 hitPixel, out vec3 hitPos, out float iterations) {
    float rayLength = ((csOrig.z + csDir.z * maxDistance) > -nearPlane) ? (-nearPlane - csOrig.z) / csDir.z : maxDistance;
    vec3 csEndPos = csOrig + csDir * rayLength;

    vec4 H0 = projToPixel * vec4(csOrig, 1.0);
    vec4 H1 = projToPixel * vec4(csEndPos, 1.0);
    float k0 = 1.0 / H0.w, k1 = 1.0 / H1.w;

    vec3 Q0 = csOrig * k0, Q1 = csEndPos * k1;
    vec2 P0 = H0.xy * k0, P1 = H1.xy * k1;
    P1 += vec2((distanceSquared(P0, P1) < 0.0001) ? 0.01 : 0.0);
    vec2 delta = P1 - P0;

    bool permute = false;
    if (abs(delta.x) < abs(delta.y)) {
        permute = true;
        delta = delta.yx;
        P0 = P0.yx;
        P1 = P1.yx;
    }

    float stepDir = sign(delta.x);
    float invdx = stepDir / delta.x;

    vec3 dQ = (Q1 - Q0) * invdx;
    float dk = (k1 - k0) * invdx;
    vec2 dP = vec2(stepDir, delta.y * invdx);

    float strideScaler = 1.0 - min(1.0, -csOrig.z / strideZCutoff);
    float pixelStride = 1.0 + strideScaler * stride;

    dP *= pixelStride; dQ *= pixelStride; dk *= pixelStride;
    P0 += dP * jitter; Q0 += dQ * jitter; k0 += dk * jitter;

    float end = P1.x * stepDir;
    float i, zA = csOrig.z, zB = csOrig.z;
    vec4 pqk = vec4(P0, Q0.z, k0);
    vec4 dPQK = vec4(dP, dQ.z, dk);
    bool intersect = false;
    for (i = 0; i < maxSteps && intersect == false && pqk.x * stepDir <= end; i++) {
        pqk += dPQK;

        zA = zB;
        zB = (dPQK.z * 0.5 + pqk.z) / (dPQK.w * 0.5 + pqk.w);

        hitPixel = permute ? pqk.yx : pqk.xy;
        hitPixel = hitPixel / screenSize;
        float currZ = linear01Depth(texture(DepthMap, hitPixel).x) * -farPlane;
        intersect = zA >= currZ - zThickness && zB <= currZ;
    }

    float addDQ = 0.0;
    if (pixelStride > 1.0 && intersect) {
        pqk -= dPQK;
        dPQK /= pixelStride;
        float originalStride = pixelStride * 0.5;
        float stride = originalStride;
        zA = pqk.z / pqk.w;
        zB = zA;
        for (float j = 0; j < binarySearchIterations; j++) {
            pqk += dPQK * stride;
            addDQ += stride;

            zA = zB;
            zB = (dPQK.z * 0.5 + pqk.z) / (dPQK.w * 0.5 + pqk.w);

            hitPixel = permute ? pqk.yx : pqk.xy;
            hitPixel = hitPixel / screenSize;
            float currZ = linear01Depth(texture(DepthMap, hitPixel).x) * -farPlane;
            bool intersect2 = zA >= currZ - zThickness && zB <= currZ;

            originalStride *= 0.5;
            stride = intersect2 ? -originalStride : originalStride;
        }
    }

    Q0.xy += dQ.xy * (i - 1) + (dQ.xy / pixelStride) * addDQ;
    Q0.z = pqk.z;
    hitPos = Q0 / pqk.w;
    iterations = i;
    return intersect;
}

void main() {
    vec2 texCoords = getTexCoords();
    float metallic = texture(AORoughMetallicMap, texCoords).b;
    if (metallic < 0.01) {
        FragColor.rgb = texture(RawFinalImage, texCoords).rgb;
        return;
    }

    float depth = texture(DepthMap, texCoords).x;
    if (depth >= 0.9999f) {
        FragColor.rgb = texture(RawFinalImage, texCoords).rgb;
        return;
    }

    vec3 normal_wS = normalize(texture(NormalMap, texCoords).xyz);
    vec3 normal_vS = mat3(sceneView) * normal_wS;

    vec4 clipPos = vec4(texCoords.xy * 2.0 - 1.0, depth * 2.0 - 1.0, 1);
    vec4 viewPos = sceneInvProjection * clipPos;
    viewPos /= viewPos.w;
    vec3 pos_vS = viewPos.xyz;
    vec3 rayDir_vS = normalize(pos_vS);
    vec3 reflectDir_vS = reflect(rayDir_vS, normal_vS);

    vec2 hitPixel = vec2(0);
    vec3 hitPos = vec3(0);
    vec2 uv = texCoords * screenSize;
    float jitter = mod((uv.x + uv.y) * 0.25, 1.0);
    float iterations = 0;
    bool hit = raymarch(pos_vS, reflectDir_vS, jitter * jitterAmount, hitPixel, hitPos, iterations);

    FragColor.rgb = texture(RawFinalImage, hitPixel).rgb;
}
