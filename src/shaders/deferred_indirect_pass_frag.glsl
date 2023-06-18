#version 420 core
out vec4 FragColor;

// SphericalGaussian(dir) := Amplitude * exp(Sharpness * (dot(Axis, dir) - 1.0f))
struct SG {
    vec3 amplitude;
    vec3 axis;
    float sharpness;
};

// AnisotropicSphericalGaussian(dir) :=
//    Amplitude * exp(-SharpnessX * dot(BasisX, dir)^2 - SharpnessY * dot(BasisY, dir)^2)
struct ASG
{
    vec3 amplitude;
    vec3 basisZ;              // Direction the ASG points
    vec3 basisX;
    vec3 basisY;
    float sharpnessX;           // Scale of the X axis
    float sharpnessY;           // Scale of the Y axis
};

uniform sampler2D PositionMap;
uniform sampler2D AlbedoMap;
uniform sampler2D AORoughMetallicMap;
uniform sampler2D NormalMap;
uniform vec2 screenSize;
uniform vec3 viewPos;

uniform ivec3 probeGridDims;
uniform vec3 gridCornerCoord;
uniform vec3 gridCellSize;
uniform int SGCount;
uniform sampler1DArray SGAmplitudes;
uniform vec3 SGDirections[12];
uniform float SGSharpness;
uniform sampler2DArray distanceOctMap;

const float Pi = 3.141592654f;

vec3 approximateSGIntegral(in SG sg);
vec3 SGIrradianceFitted(in SG lightingLobe, in vec3 normal);
ASG warpDistributionASG(in SG ndf, in vec3 view);
vec3 specularTermASGWarp(in SG light, in vec3 normal, in float roughness, in vec3 view, in vec3 specAlbedo);
void computeSGContribution(in int probeIndex, in vec3 normal, in vec3 specularAlbedo, in float roughness, in vec3 view, out vec3 irradiance, out vec3 specular);

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
    vec3 diffuseAlbedo = mix(albedo, vec3(0.0), vec3(metallic));

    ivec3 baseCoord = baseGridCoord(worldPos);
    vec3 basePos = gridCoordToPosition(baseCoord);
    vec3 alpha = clamp((worldPos - basePos) / gridCellSize, vec3(0.0), vec3(1.0));  // interpolation weights

    // selects 8 probes surrounding fragment
    vec3 irradiance = vec3(0.0);
    float sum_weight = 0.0;
    for (int idx = 0; idx < 8; idx++) {
        ivec3 offset = ivec3(idx, idx >> 1, idx >> 2) & ivec3(1);
        ivec3 probeGridCoord = ivec3(clamp(baseCoord + offset, vec3(0.0), vec3(probeGridDims.x, probeGridDims.y, probeGridDims.z) - vec3(1.0)));
        int probeIndex = gridCoordToFlatIndex(vec3(probeGridCoord));

        vec3 trilinear = mix(1.0 - alpha, alpha, vec3(offset));
        float weight = trilinear.x * trilinear.y * trilinear.z;

        vec3 probePos = gridCoordToPosition(probeGridCoord);
        vec3 probeToPos = worldPos - probePos;
        vec3 dir = normalize(-probeToPos);

        weight *= max(0.05, dot(dir, normal));

        vec2 octDir = octEncode(-dir) * 0.5 + 0.5;
        vec2 temp = texture(distanceOctMap, vec3(octDir, probeIndex)).xy;
        float mean = temp.x;
        float variance = abs(temp.y - (mean * mean));

        float probeDist = length(probeToPos);
        float t_sub_mean = probeDist - mean;
        float chebychev = variance / (variance + (t_sub_mean * t_sub_mean));
        weight *= ((probeDist <= mean)) ? 1.0 : max(chebychev, 0.0);
        irradiance = vec3(mean, variance, probeDist);

        weight = max(0.002, weight);
        sum_weight += weight;

        vec3 indirectDiffuse = vec3(0.0);
        vec3 indirectSpecular = vec3(0.0);
        computeSGContribution(probeIndex, normal, F0, roughness, viewDir, indirectDiffuse, indirectSpecular);

        irradiance += indirectDiffuse * (diffuseAlbedo / Pi) * weight;
    }

    FragColor.rgb = irradiance * ao / sum_weight;
}

//-------------------------------------------------------------------------------------------------
// Computes the inner product of two SG's, which is equal to Integrate(SGx(v) * SGy(v) * dv).
//-------------------------------------------------------------------------------------------------
vec3 SGInnerProduct(in SG x, in SG y) {
    float umLength = length(x.sharpness * x.axis + y.sharpness * y.axis);
    vec3 expo = exp(umLength - x.sharpness - y.sharpness) * x.amplitude * y.amplitude;
    float other = 1.0f - exp(-2.0f * umLength);
    return (2.0f * Pi * expo * other) / umLength;
}

//-------------------------------------------------------------------------------------------------
// Computes the approximate integral of an SG over the entire sphere. The error vs. the
// non-approximate version decreases as sharpeness increases.
//-------------------------------------------------------------------------------------------------
vec3 approximateSGIntegral(in SG sg) {
    return 2 * Pi * (sg.amplitude / sg.sharpness);
}

//-------------------------------------------------------------------------------------------------
// Computes the approximate incident irradiance from a single SG lobe containing incoming radiance.
// The irradiance is computed using a fitted approximation polynomial. This approximation
// and its implementation were provided by Stephen Hill.
//-------------------------------------------------------------------------------------------------
vec3 SGIrradianceFitted(in SG lightingLobe, in vec3 normal) {
    const float muDotN = dot(lightingLobe.axis, normal);
    const float lambda = lightingLobe.sharpness;

    const float c0 = 0.36;
    const float c1 = 1.0 / (4.0 * c0);

    float eml  = exp(-lambda);
    float em2l = eml * eml;
    float rl   = 1.0 / lambda;

    float scale = 1.0 + 2.0 * em2l - rl;
    float bias  = (eml - em2l) * rl - em2l;

    float x  = sqrt(1.0 - scale);
    float x0 = c0 * muDotN;
    float x1 = c1 * x;

    float n = x0 + x1;

    float y = (abs(x0) <= x1) ? n * n / x : clamp(muDotN, 0, 1);

    float normalizedIrradiance = scale * y + bias;

    return normalizedIrradiance * approximateSGIntegral(lightingLobe);
}

//-------------------------------------------------------------------------------------------------
// Returns an SG approximation of the GGX NDF used in the specular BRDF. For a single-lobe
// approximation, the resulting NDF actually more closely resembles a Beckmann NDF.
//-------------------------------------------------------------------------------------------------
SG distributionTermSG(in vec3 direction, in float roughness) {
    SG distribution;
    distribution.axis = direction;
    float m2 = roughness * roughness;
    distribution.sharpness = 2 / m2;
    distribution.amplitude = vec3(1.0 / (Pi * m2));

    return distribution;
}

//-------------------------------------------------------------------------------------------------
// Generate an ASG that best represents the NDF SG but with it's axis oriented in the direction
// of the current BRDF slice. This will allow easier integration, because the SG\ASG are both
// in the same domain.
//
// The warped NDF can be represented better as an ASG, so following Kun Xu from
// 'Anisotropic Spherical Gaussians' we change the SG to an ASG because the distribution of
// an NDF stretches at grazing angles.
//-------------------------------------------------------------------------------------------------
ASG warpDistributionASG(in SG ndf, in vec3 view) {
    ASG warp;

    // Generate any orthonormal basis with Z pointing in the direction of the reflected view vector
    warp.basisZ = reflect(-view, ndf.axis);
    warp.basisX = normalize(cross(ndf.axis, warp.basisZ));
    warp.basisY = normalize(cross(warp.basisZ, warp.basisX));

    warp.amplitude = ndf.amplitude;

    float VdotA = max(dot(view, ndf.axis), 0.1);
    warp.sharpnessX = ndf.sharpness / (8.0 * VdotA * VdotA);
    warp.sharpnessY = ndf.sharpness / 8.0;

    return warp;
}

float GGX_V1(float NdotV, float m2) {
    return 1.0 / (NdotV + sqrt(m2 + (1 - m2) * NdotV * NdotV));
}

//-------------------------------------------------------------------------------------------------
// Evaluates an ASG given a direction on a unit sphere
//-------------------------------------------------------------------------------------------------
vec3 evaluateASG(in ASG asg, in vec3 dir) {
    float smoothTerm = clamp(dot(asg.basisZ, dir), 0, 1);
    float lambdaTerm = asg.sharpnessX * dot(dir, asg.basisX) * dot(dir, asg.basisX);
    float muTerm = asg.sharpnessY * dot(dir, asg.basisY) * dot(dir, asg.basisY);
    return asg.amplitude * smoothTerm * exp(-lambdaTerm - muTerm);
}

//-------------------------------------------------------------------------------------------------
// Convolve an SG with an ASG
//-------------------------------------------------------------------------------------------------
vec3 convolveASG_SG(in ASG asg, in SG sg) {
    // The ASG paper specifes an isotropic SG as exp(2 * nu * (dot(v, axis) - 1)),
    // so we must divide our SG sharpness by 2 in order to get the nup parameter expected by
    // the ASG formulas
    float nu = sg.sharpness * 0.5f;

    ASG convolveASG;
    convolveASG.basisX = asg.basisX;
    convolveASG.basisY = asg.basisY;
    convolveASG.basisZ = asg.basisZ;

    convolveASG.sharpnessX = (nu * asg.sharpnessX) / (nu + asg.sharpnessX);
    convolveASG.sharpnessY = (nu * asg.sharpnessY) / (nu + asg.sharpnessY);

    convolveASG.amplitude = vec3(Pi / sqrt((nu + asg.sharpnessX) * (nu + asg.sharpnessY)));

    return evaluateASG(convolveASG, sg.axis) * sg.amplitude * asg.amplitude;
}

//-------------------------------------------------------------------------------------------------
// Computes the specular reflectance from a single SG lobe containing incoming radiance
//-------------------------------------------------------------------------------------------------
vec3 specularTermASGWarp(in SG light, in vec3 normal, in float roughness, in vec3 view, in vec3 specAlbedo) {
    // Create an SG that approximates the NDF. Note that a single SG lobe is a poor fit for
    // the GGX NDF, since the GGX distribution has a longer tail. A sum of 3 SG's can more
    // closely match the shape of a GGX distribution, but it would also increase the cost
    // computing specular by a factor of 3.
    SG ndf = distributionTermSG(normal, roughness);

    // Apply a warpring operation that will bring the SG from the half-angle domain the the
    // the lighting domain. The resulting lobe is another SG.
    ASG warpedNDF = warpDistributionASG(ndf, view);

    // Convolve the NDF with the SG light
    vec3 result = convolveASG_SG(warpedNDF, light);

    // Parameters needed for evaluating the visibility term
    float m2 = roughness * roughness;
    float NdotL = clamp(dot(normal, warpedNDF.basisZ), 0, 1);
    float NdotV = clamp(dot(normal, view), 0, 1);
    vec3 h = normalize(warpedNDF.basisZ + view);

    // The visibility term is evaluated at the center of our warped BRDF lobe
    result *= GGX_V1(m2, NdotL) * GGX_V1(m2, NdotV);

    // Fresnel evaluated at the center of our warped BRDF lobe
    result *= specAlbedo + (1.0 - specAlbedo) * pow((1.0 - clamp(dot(warpedNDF.basisZ, h), 0, 1)), 5.0);

    // Fade out spec entirely when lower than 0.1% albedo
    result *= clamp(dot(specAlbedo, vec3(333.0)), 0, 1);

    // Cosine term evaluated at the center of our warped BRDF lobe
    result *= NdotL;

    return max(result, 0.0);
}

// ------------------------------------------------------------------------------------------------
// Determine the exit radiance towards the eye from the SG's stored in the lightmap
// ------------------------------------------------------------------------------------------------
void computeSGContribution(in int probeIndex, in vec3 normal, in vec3 specularAlbedo,
in float roughness, in vec3 view, out vec3 irradiance, out vec3 specular) {
    irradiance = vec3(0.0);
    specular = vec3(0.0);

    for (int i = 0; i < SGCount; i++) {
        SG sg;
        sg.amplitude = texture(SGAmplitudes, vec2(i, probeIndex)).xyz;
        sg.axis = SGDirections[i].xyz;
        sg.sharpness = SGSharpness;

        irradiance += SGIrradianceFitted(sg, normal);
        specular += specularTermASGWarp(sg, normal, roughness, view, specularAlbedo);
    }
}
