#version 420 core
out vec4 FragColor;

uniform sampler2D DepthMap;
uniform sampler2D NormalMap;
uniform sampler2D AORoughMetallicMap;
uniform sampler2D RawFinalImage;

uniform vec2 screenSize;
uniform mat4 sceneView;
uniform mat4 sceneProjection;
uniform mat4 sceneInvProjection;

const float rayStep = 0.1;
const float maxSteps = 100;
const float distanceBias = 0.05;

vec2 getTexCoords() {
    return gl_FragCoord.xy / screenSize;
}

vec3 positionFromDepth(float depth, vec2 texCoords) {
    vec4 pos_sS = vec4(texCoords * 2.0 - 1.0, depth, 1.0);
    vec4 pos_vS = sceneInvProjection * pos_sS;
    pos_vS /= pos_vS.w;

    return pos_vS.xyz;
}

vec2 projectedPosition(vec3 pos) {
    vec4 samplePosition = sceneProjection * vec4(pos, 1.f);
    samplePosition.xy = (samplePosition.xy / samplePosition.w) * 0.5 + 0.5;
    return samplePosition.xy;
}

vec3 SSR(vec3 position, vec3 reflection) {
    vec3 stepSize = rayStep * reflection;
    vec3 marchingPos = position + stepSize;
    float delta;
    float depthFromScreen;
    vec2 screenPos;

    int i = 0;
    for (; i < maxSteps; i++) {
        screenPos = projectedPosition(marchingPos);
        depthFromScreen = abs(positionFromDepth(texture(DepthMap, screenPos).x, screenPos).z);
        delta = abs(marchingPos.z) - depthFromScreen;

        if (abs(delta) < distanceBias) {
            return texture(RawFinalImage, screenPos).xyz;
        }
        if (delta > 0) {
            break;
        }
        marchingPos += stepSize;
    }
    for (; i < maxSteps; i++) {
        stepSize *= 0.5;
        marchingPos = marchingPos - stepSize * sign(delta);

        screenPos = projectedPosition(marchingPos);
        depthFromScreen = abs(positionFromDepth(texture(DepthMap, screenPos).x, screenPos).z);
        delta = abs(marchingPos.z) - depthFromScreen;

        if (abs(delta) < distanceBias) {
            return texture(RawFinalImage, screenPos).xyz;
        }
    }

    return vec3(0.0);
}

void main() {
    vec2 texCoords = getTexCoords();
    float metallic = texture(AORoughMetallicMap, texCoords).b;
    if (metallic < 0.01) {
        discard;
    }

    float depth = texture(DepthMap, texCoords).x;
    vec3 normal = texture(NormalMap, texCoords).xyz;

    // position + normal in view space
    vec3 pos_vS = positionFromDepth(depth, texCoords);
    vec3 normal_vS = (sceneView * vec4(normal, 0.0)).xyz;

    vec3 reflectedDir = normalize(reflect(pos_vS, normalize(normal_vS)));
    vec3 color = SSR(pos_vS, reflectedDir);

    FragColor.rgb = color;
}
