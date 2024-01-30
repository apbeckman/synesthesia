vec3 iResolution = vec3(RENDERSIZE.x, RENDERSIZE.y, 0.0);
vec3 iChannelResolution = iResolution;
float iGlobalTime = TIME*1.0;
int iFrame = int(FRAMECOUNT);
float time = TIME;
vec2 resolution = RENDERSIZE;
vec4 iMouse = vec4(0.5);
float dist(vec2 p0, vec2 pf){return sqrt((pf.x-p0.x)*(pf.x-p0.x)+(pf.y-p0.y)*(pf.y-p0.y));}
float d = dist(RENDERSIZE.xy*0.5,_xy.xy)*(_mouse.x/RENDERSIZE.x+0.1)*0.005;
float d2 = dist(RENDERSIZE.xy*0.5,_xy.xy)*0.0125;


// ****************** PASS 0 ***********************
#define pi2_inv 0.159154943091895335768883763372

vec2 lower_left(vec2 uv)
{
    return fract(uv * 0.5);
}

vec2 lower_right(vec2 uv)
{
    return fract((uv - vec2(1, 0.)) * 0.5);
}

vec2 upper_left(vec2 uv)
{
    return fract((uv - vec2(0., 1)) * 0.5);
}

vec2 upper_right(vec2 uv)
{
    return fract((uv - 1.) * 0.5);
}

vec4 BlurA(vec2 uv, int level)
{
    if(level <= 0)
    {
        return texture(buffA, fract(uv));
    }

    uv = upper_left(uv);
    for(int depth = 1; depth < 8; depth++)
    {
        if(depth >= level)
        {
            break;
        }
        uv = lower_right(uv);
    }

    return texture(buffD, uv);
}

vec2 GradientA(vec2 uv, vec2 d, vec4 selector, int level){
    vec4 dX = 0.5*BlurA(uv + vec2(1.,0.)*d, level) - 0.5*BlurA(uv - vec2(1.,0.)*d, level);
    vec4 dY = 0.5*BlurA(uv + vec2(0.,1.)*d, level) - 0.5*BlurA(uv - vec2(0.,1.)*d, level);
    dY += shock;
    return vec2( dot(dX, selector), dot(dY, selector) );
}

vec2 rot90(vec2 vector){
    return vector.yx*vec2(1,-1);
}

vec2 complex_mul(vec2 factorA, vec2 factorB){
    return vec2( factorA.x*factorB.x - factorA.y*factorB.y, factorA.x*factorB.y + factorA.y*factorB.x);
}

float sigmoid(float x) {
    return 2./(1. + exp2(-x)) - 1.;
}

float conetip(vec2 uv, vec2 pos, float size, float min)
{
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);
    return max( min, 1. - length((uv - pos) * aspect / size) );
}

float warpFilter(vec2 uv, vec2 pos, float size, float ramp)
{
    return 0.5 + sigmoid( conetip(uv, pos, size, -16.) * ramp) * 0.5;
}

vec2 vortex_warp(vec2 uv, vec2 pos, float size, float ramp, vec2 rot)
{
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);

    vec2 pos_correct = 0.5 + (pos - 0.5);
    vec2 rot_uv = pos_correct + complex_mul((uv - pos_correct)*aspect, rot)/aspect;
    float f1lter = warpFilter(uv, pos_correct, size, ramp);
    return mix(uv, rot_uv, f1lter);
}

vec2 vortex_pair_warp(vec2 uv, vec2 pos, vec2 vel)
{
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);
    float ramp = 4.;

    float d = 0.425;

    float l = length(vel);
    vec2 p1 = pos;
    vec2 p2 = pos;

    if(l > 0.){
        vec2 normal = normalize(vel.yx * vec2(-1., 1.))/aspect;
        p1 = pos - normal * d / 2.;
        p2 = pos + normal * d / 2.;
    }

    float w = l / d * 2.;

    // two overlapping rotations that would annihilate when they were not displaced.
    vec2 circle1 = vortex_warp(uv, p1, d, ramp, vec2(cos(w),sin(w)));
    vec2 circle2 = vortex_warp(uv, p2, d, ramp, vec2(cos(-w),sin(-w)));
    return (circle1 + circle2) / 2.;
}

void mainImage0( out vec4 fragColor, in vec2 fragCoord )
{
    fragCoord -= _uvc*PI*Zoom;
    fragCoord += _uvc*Stretch;
    fragCoord += _uvc*fisheye*d2;
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec4 noise = (texture(colornoise, _uv*4.0 + fract(vec2(42,56)*TIME*0.125))-0.5)*2.;
    uv += _uvc*pow(tunnel,vec2(4.0));

    if ((pyramid>0.5)&&(uv.x > uv.y*0.5)){
        uv.y = 0.5-uv.y;
    }


    if(iFrame<10)
    {
        fragColor = noise;
        return;
    }

    float intensity = pow(syn_Presence*0.5+syn_Intensity*0.5,5.0);
    vec2 flowVel = vec2(0.00,mix(-0.0001, -0.002, intensity));
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);
    vec2 pixelSize = 1. / iResolution.xy;

    uv = vortex_pair_warp(uv, vec2(0.5, 0.5), flowVel*aspect*1.4);

    // expansion
    vec2 gradientLookupDistance = pixelSize*5.;
    float expansionFactor = mix(0.3, 1.5, syn_Intensity*0.5+syn_Level*0.5);
    // bool jumpyExpoRegime = false;
    // if (jumpyExpoRegime == true){
    //     expansionFactor = 1.*syn_OnBeat*4.0;
    // }


    // reaction-diffusion
    float differentialFactor = 12./256.;
    float increment = -3./256.;
    float noiseFactor = 2./256.;

    // rock-paper-scissor
    float feedBack = shimmer*6./256.;
    float feedForward = shimmer*6./256.;

    float motion = clamp(-high+bass,-0.35,1.0)*0.01;
    if (bpm_oscillate > 0.5){
        // expansionFactor = syn_BPMSin2;
        motion = clamp(syn_BPMSin2*2-1.0,-1.0,1.0)*0.01+0.002*syn_HighHits;
    }

    uv += GradientA(uv, gradientLookupDistance, vec4(0., 0., 0., 1.), 1)*motion;

    fragColor.r = BlurA(uv + GradientA(uv, gradientLookupDistance, vec4(4.,0.,-2.,0.), 1)*pixelSize*expansionFactor, 0).r;
    fragColor.g = BlurA(uv + GradientA(uv, gradientLookupDistance, vec4(0.,4.,0.,-2.), 1)*pixelSize*expansionFactor, 0).g;
    fragColor.b = BlurA(uv + GradientA(uv, gradientLookupDistance, vec4(-2.,0.,4.,0.), 1)*pixelSize*expansionFactor, 0).b;
    fragColor.a = BlurA(uv + GradientA(uv, gradientLookupDistance, vec4(0.,-2.,0.,4.), 1)*pixelSize*expansionFactor, 0).a;

    if (_exists(syn_UserImage)){
        float userImg = _loadUserImageAsMask().r;
        if (media_color<1.0){
            fragColor.r += userImg;
        } else if (media_color<2.0){
            fragColor.g += userImg;
        } else if (media_color<3.0){
            fragColor.b += userImg;
        } else if (media_color<4.0){
            fragColor.a += userImg;
        }
    }

    fragColor += (BlurA(uv, 1) - BlurA(uv, 2))*differentialFactor;

    fragColor += increment + noise * noiseFactor;

    fragColor = fragColor + (1.0-smoothstep(0.1,0.4,length(_uvc)))*reset_center*(0.8+0.2*noise*0.5);

    if (length(_uv)<0.01){
        fragColor = vec4(1.0);
    }

    fragColor -= fragColor.argb * feedBack;
    fragColor += fragColor.gbar * feedForward;

    fragColor = clamp(fragColor, 0., 1.);

//    fragColor = noise; // reset
}

// *************** PASS 2 *****************
vec4 blur_horizontal(sampler2D channel, vec2 uv, float scale)
{
    float h = scale / iResolution.x;
    vec4 sum = vec4(0.0);

    sum += texture(channel, fract(vec2(uv.x - 4.0*h, uv.y)) ) * 0.05;
    sum += texture(channel, fract(vec2(uv.x - 3.0*h, uv.y)) ) * 0.09;
    sum += texture(channel, fract(vec2(uv.x - 2.0*h, uv.y)) ) * 0.12;
    sum += texture(channel, fract(vec2(uv.x - 1.0*h, uv.y)) ) * 0.15;
    sum += texture(channel, fract(vec2(uv.x + 0.0*h, uv.y)) ) * 0.16;
    sum += texture(channel, fract(vec2(uv.x + 1.0*h, uv.y)) ) * 0.15;
    sum += texture(channel, fract(vec2(uv.x + 2.0*h, uv.y)) ) * 0.12;
    sum += texture(channel, fract(vec2(uv.x + 3.0*h, uv.y)) ) * 0.09;
    sum += texture(channel, fract(vec2(uv.x + 4.0*h, uv.y)) ) * 0.05;

    return sum/0.98; // normalize
}

vec4 blur_horizontal_left_column(vec2 uv, int depth)
{
    float h = pow(2., float(depth)) / iResolution.x;
    vec2 uv1, uv2, uv3, uv4, uv5, uv6, uv7, uv8, uv9;

    uv1 = fract(vec2(uv.x - 4.0 * h, uv.y) * 2.);
    uv2 = fract(vec2(uv.x - 3.0 * h, uv.y) * 2.);
    uv3 = fract(vec2(uv.x - 2.0 * h, uv.y) * 2.);
    uv4 = fract(vec2(uv.x - 1.0 * h, uv.y) * 2.);
    uv5 = fract(vec2(uv.x + 0.0 * h, uv.y) * 2.);
    uv6 = fract(vec2(uv.x + 1.0 * h, uv.y) * 2.);
    uv7 = fract(vec2(uv.x + 2.0 * h, uv.y) * 2.);
    uv8 = fract(vec2(uv.x + 3.0 * h, uv.y) * 2.);
    uv9 = fract(vec2(uv.x + 4.0 * h, uv.y) * 2.);

    if(uv.y > 0.5)
    {
        uv1 = upper_left(uv1);
        uv2 = upper_left(uv2);
        uv3 = upper_left(uv3);
        uv4 = upper_left(uv4);
        uv5 = upper_left(uv5);
        uv6 = upper_left(uv6);
        uv7 = upper_left(uv7);
        uv8 = upper_left(uv8);
        uv9 = upper_left(uv9);
    }
    else{
        uv1 = lower_left(uv1);
        uv2 = lower_left(uv2);
        uv3 = lower_left(uv3);
        uv4 = lower_left(uv4);
        uv5 = lower_left(uv5);
        uv6 = lower_left(uv6);
        uv7 = lower_left(uv7);
        uv8 = lower_left(uv8);
        uv9 = lower_left(uv9);
    }

    for(int level = 0; level < 8; level++)
    {
        if(level >= depth)
        {
            break;
        }

        uv1 = lower_right(uv1);
        uv2 = lower_right(uv2);
        uv3 = lower_right(uv3);
        uv4 = lower_right(uv4);
        uv5 = lower_right(uv5);
        uv6 = lower_right(uv6);
        uv7 = lower_right(uv7);
        uv8 = lower_right(uv8);
        uv9 = lower_right(uv9);
    }

    vec4 sum = vec4(0.0);

    sum += texture(buffD, uv1) * 0.05;
    sum += texture(buffD, uv2) * 0.09;
    sum += texture(buffD, uv3) * 0.12;
    sum += texture(buffD, uv4) * 0.15;
    sum += texture(buffD, uv5) * 0.16;
    sum += texture(buffD, uv6) * 0.15;
    sum += texture(buffD, uv7) * 0.12;
    sum += texture(buffD, uv8) * 0.09;
    sum += texture(buffD, uv9) * 0.05;

    return sum/0.98; // normalize
}

void mainImage2( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / iResolution.xy;

    if(uv.x < 0.5)
    {
        vec2 uv_half = fract(uv*2.);
        if(uv.y > 0.5)
        {
            fragColor = blur_horizontal(buffA, uv_half, 1.);
        }
        else
        {
            fragColor = blur_horizontal(buffA, uv_half, 1.);
        }
    }
    else
    {
        for(int level = 0; level < 8; level++)
        {
            if((uv.x > 0.5 && uv.y > 0.5) || (uv.x <= 0.5))
            {
                break;
            }
            vec2 uv_half = fract(uv*2.);
            fragColor = blur_horizontal_left_column(uv_half, level);
            uv = uv_half;
        }
    }
}


// *************** PASS 4 *****************

vec4 blur_vertical_upper_left(sampler2D channel, vec2 uv)
{
    float v = 1. / iResolution.y;
    vec4 sum = vec4(0.0);
    sum += texture(channel, upper_left(vec2(uv.x, uv.y - 4.0*v)) ) * 0.05;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y - 3.0*v)) ) * 0.09;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y - 2.0*v)) ) * 0.12;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y - 1.0*v)) ) * 0.15;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 0.0*v)) ) * 0.16;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 1.0*v)) ) * 0.15;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 2.0*v)) ) * 0.12;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 3.0*v)) ) * 0.09;
    sum += texture(channel, upper_left(vec2(uv.x, uv.y + 4.0*v)) ) * 0.05;
    return sum/0.98; // normalize
}

vec4 blur_vertical_lower_left(sampler2D channel, vec2 uv)
{
    float v = 1. / iResolution.y;
    vec4 sum = vec4(0.0);
    sum += texture(channel, lower_left(vec2(uv.x, uv.y - 4.0*v)) ) * 0.05;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y - 3.0*v)) ) * 0.09;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y - 2.0*v)) ) * 0.12;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y - 1.0*v)) ) * 0.15;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 0.0*v)) ) * 0.16;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 1.0*v)) ) * 0.15;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 2.0*v)) ) * 0.12;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 3.0*v)) ) * 0.09;
    sum += texture(channel, lower_left(vec2(uv.x, uv.y + 4.0*v)) ) * 0.05;
    return sum/0.98; // normalize
}

vec4 blur_vertical_left_column(vec2 uv, int depth)
{
    float v = pow(2., float(depth)) / iResolution.y;

    vec2 uv1, uv2, uv3, uv4, uv5, uv6, uv7, uv8, uv9;

    uv1 = fract(vec2(uv.x, uv.y - 4.0*v) * 2.);
    uv2 = fract(vec2(uv.x, uv.y - 3.0*v) * 2.);
    uv3 = fract(vec2(uv.x, uv.y - 2.0*v) * 2.);
    uv4 = fract(vec2(uv.x, uv.y - 1.0*v) * 2.);
    uv5 = fract(vec2(uv.x, uv.y + 0.0*v) * 2.);
    uv6 = fract(vec2(uv.x, uv.y + 1.0*v) * 2.);
    uv7 = fract(vec2(uv.x, uv.y + 2.0*v) * 2.);
    uv8 = fract(vec2(uv.x, uv.y + 3.0*v) * 2.);
    uv9 = fract(vec2(uv.x, uv.y + 4.0*v) * 2.);

    if(uv.y > 0.5)
    {
        uv1 = upper_left(uv1);
        uv2 = upper_left(uv2);
        uv3 = upper_left(uv3);
        uv4 = upper_left(uv4);
        uv5 = upper_left(uv5);
        uv6 = upper_left(uv6);
        uv7 = upper_left(uv7);
        uv8 = upper_left(uv8);
        uv9 = upper_left(uv9);
    }
    else{
        uv1 = lower_left(uv1);
        uv2 = lower_left(uv2);
        uv3 = lower_left(uv3);
        uv4 = lower_left(uv4);
        uv5 = lower_left(uv5);
        uv6 = lower_left(uv6);
        uv7 = lower_left(uv7);
        uv8 = lower_left(uv8);
        uv9 = lower_left(uv9);
    }

    for(int level = 0; level < 8; level++)
    {
        if(level > depth)
        {
            break;
        }

        uv1 = lower_right(uv1);
        uv2 = lower_right(uv2);
        uv3 = lower_right(uv3);
        uv4 = lower_right(uv4);
        uv5 = lower_right(uv5);
        uv6 = lower_right(uv6);
        uv7 = lower_right(uv7);
        uv8 = lower_right(uv8);
        uv9 = lower_right(uv9);
    }

    vec4 sum = vec4(0.0);

    sum += texture(buffC, uv1) * 0.05;
    sum += texture(buffC, uv2) * 0.09;
    sum += texture(buffC, uv3) * 0.12;
    sum += texture(buffC, uv4) * 0.15;
    sum += texture(buffC, uv5) * 0.16;
    sum += texture(buffC, uv6) * 0.15;
    sum += texture(buffC, uv7) * 0.12;
    sum += texture(buffC, uv8) * 0.09;
    sum += texture(buffC, uv9) * 0.05;

    return sum/0.98; // normalize
}

void mainImage3( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec2 uv_orig = uv;
    vec2 uv_half = fract(uv*2.);
    if(uv.x < 0.5)
    {
        if(uv.y > 0.5)
        {
            fragColor = blur_vertical_upper_left(buffC, uv_half);
        }
        else
        {
            fragColor = blur_vertical_lower_left(buffC, uv_half);
        }
    }
    else
    {
        for(int level = 0; level < 8; level++)
        {
            if((uv.x > 0.5 && uv.y >= 0.5) || (uv.x < 0.5))
            {
                break;
            }
            vec2 uv_half = fract(uv*2.);
            fragColor = blur_vertical_left_column(uv_half, level);
            uv = uv_half;
        }
    }
    uv = uv_orig;
    float eighth = 1./8.;
    if(uv.x > 7.*eighth && uv.x < 8.*eighth && uv.y > 2.*eighth && uv.y < 3.*eighth)
    {
        fragColor = vec4(iMouse.xy / iResolution.xy, iMouse.zw / iResolution.xy);
    }
}

// *************** Main Image Pass *****************

void mainImage4( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = _uv;

    // if (uv.x*mod(syn_BeatTime*5.0+syn_RandomOnBeat*2.0,8.0)*0.5+uv.y*mod(syn_BeatTime*13.0+syn_RandomOnBeat*1.7,8.0)*0.5 > 1.0){
    //     uv = 1.0-uv;
    // }
    if ((uv.x > 0.5)&&(mirror>0.5)){
        uv.x = 1.0-uv.x;
    }

    vec2 pixelSize = 1. / iResolution.xy;
    vec2 aspect = vec2(1.,iResolution.y/iResolution.x);

    vec2 d = pixelSize*2.;
    vec4 dx = (BlurA(uv + vec2(1,0)*d, 1) - BlurA(uv - vec2(1,0)*d, 1))*0.5;
    vec4 dy = (BlurA(uv + vec2(0,1)*d, 1) - BlurA(uv - vec2(0,1)*d, 1))*0.5;

    d = pixelSize*1.;
    dx += BlurA(uv + vec2(1,0)*d, 0) - BlurA(uv - vec2(1,0)*d, 0);
    dy += BlurA(uv + vec2(0,1)*d, 0) - BlurA(uv - vec2(0,1)*d, 0);
    vec2 lightSize=vec2(0.5);

    //for B
    //midHighCol is general color?
    //bassCol1 is a color in between coral
    //highCol1 does almost nothing. creeps in on the edges.
    //bassCol2 is general color "on" the life spots
    //highCol2 does almost nothing when sim is running. creeps in from edges.
    //midCol is squiggles layer... Its red and white tho.

    //for A
    //midHighCol is sheen near the wavefronts. like "newSpawn"
    //bassCol1 is big wavefronts and smooth shading
    //highCol1 is wavefronts
    //bassCol2 crisp wavefronts and smooth wavefronts, of a diff color probably.
    //highCol2 sort of ghostly messy wavefronts
    //midCol sort of a networky layer. Looks brainular.


    vec3 midHighCol = vec3(0.7,1.66,2.0);
    vec3 bassCol1 = vec3(8.0,6.,2.);
    vec3 highCol1 = vec3(0.1,0.,0.4);
    vec3 bassCol2 = vec3(1.25,1.35,1.4);
    vec3 highCol2 = vec3(0.25,0.75,1.);
    vec3 midCol = vec3(1.,1.25,1.5);


    // float cReg1 = 0.0;
    // float cReg2 = 0.0;
    // float cReg3 = 0.0;
    // float cReg4 = 0.0;


    // *** Color Regime 1 ***
    midHighCol = mix(midHighCol, _normalizeRGB(252, 255, 245)*2.0, cReg1);
    bassCol1 = mix(bassCol1, _normalizeRGB(0, 38, 28)*2.0, cReg1);
    highCol1 = mix(highCol1, _normalizeRGB(4, 76, 41)*2.0, cReg1);
    bassCol2 = mix(bassCol2, _normalizeRGB(22, 127, 57)*2.0, cReg1);
    highCol2 = mix(highCol2, _normalizeRGB(150, 237, 137)*2.0, cReg1);
    midCol = mix(midCol, _normalizeRGB(4, 76, 41)*2.0, cReg1);
    // highCol2 = mix(highCol2, vec3(0.3,0.4,0.22), cReg1);

    // *** Color Regime 2 ***
    midHighCol = mix(midHighCol, _normalizeRGB(70, 137, 102)*2.0, cReg2);
    bassCol1 = mix(bassCol1, _normalizeRGB(255, 240, 165)*2.0, cReg2);
    highCol1 = mix(highCol1, _normalizeRGB(255, 176, 59)*2.0, cReg2);
    bassCol2 = mix(bassCol2, _normalizeRGB(182, 73, 38)*2.0, cReg2);
    highCol2 = mix(highCol2, _normalizeRGB(142, 40, 0)*2.0, cReg2);
    midCol = mix(midCol, _normalizeRGB(142, 40, 0)*2.0, cReg2);

    // *** Color Regime 3 ***
    midHighCol = mix(midHighCol, _normalizeRGB(28, 24, 34)*3.0, cReg3);
    bassCol1 = mix(bassCol1, _normalizeRGB(64, 21, 42)*3.0, cReg3);
    highCol1 = mix(highCol1, _normalizeRGB(115, 22, 48)*2.0, cReg3);
    bassCol2 = mix(bassCol2, _normalizeRGB(255, 84, 52)*2.0, cReg3);
    highCol2 = mix(highCol2, _normalizeRGB(204, 30, 44)*2.0, cReg3);
    midCol = mix(midCol, _normalizeRGB(64, 21, 42)*2.0, cReg3);

    // *** Color Regime 4 ***
    midHighCol = mix(midHighCol, _normalizeRGB(46, 9, 39)*4.0, cReg4);
    bassCol1 = mix(bassCol1, _normalizeRGB(217, 0, 0)*2.0, cReg4);
    highCol1 = mix(highCol1, _normalizeRGB(4, 117, 111)*2.0, cReg4);
    bassCol2 = mix(bassCol2, _normalizeRGB(255, 45, 0)*2.0, cReg4);
    highCol2 = mix(highCol2, _normalizeRGB(255, 140, 0)*2.0, cReg4);
    midCol = mix(midCol, _normalizeRGB(4, 117, 111)*2.0, cReg4);


    vec2 zoneMids = vec2(sin(smooth_midtime*0.23), cos(smoothTimeC*0.17));
    vec2 zoneMidHighs = vec2(sin(smooth_hightime*0.23), cos(smooth_midtime*0.17));
    vec2 zoneHighs = vec2(sin(smooth_hightime*0.23), cos(smoothTimeB*0.17));

    midHighCol *= (0.05+0.95*syn_MidHighPresence*distance(zoneMidHighs, _uvc)*0.8);
    bassCol1 *= (0.05+0.95*(sqrt(syn_BassPresence))*1.25);
    highCol1 *= (0.05+0.95*sqrt(syn_HighPresence));
    bassCol2 *= (0.05+0.95*sqrt(syn_BassPresence)*1.5);
    highCol2 *= (0.05+0.95*syn_HighHits*(0.5+syn_Presence*0.5)*distance(zoneHighs, _uvc));
    midCol *= (0.05+0.95*distance(zoneMids, _uvc)*sqrt(syn_MidPresence));

    float lowIntensity = 1.0-syn_Presence;

    vec3 finalColor = vec3(0.0);

    finalColor += BlurA(uv+vec2(dx.x,dy.x)*pixelSize*8., 0).x * midHighCol;
    finalColor = mix(finalColor, bassCol1, BlurA(uv + vec2(dx.x,dy.x)*lightSize, 3).y*0.4*0.75*vec3(1.-BlurA(uv+vec2(dx.x,dy.x)*pixelSize*8., 0).x));
    finalColor = mix(finalColor, highCol1, BlurA(uv, 1).a*length(GradientA(uv, pixelSize*2., vec4(0.,0.,0., 1.), 0))*5.);
    finalColor = mix(finalColor, bassCol2, BlurA(uv, 0).x*BlurA(uv + GradientA(uv, pixelSize*2.5, vec4(-256.,32.,-128.,32.), 1)*pixelSize, 2).y);
    finalColor = mix(finalColor, highCol2, BlurA(uv, 1).x*length(GradientA(uv+GradientA(uv, pixelSize*2., vec4(0.,0.,128.,0.), 1)*pixelSize, pixelSize*2., vec4(0.,0.,0.,1.), 0))*5.);
    finalColor = mix(finalColor, midCol, 0.5*(1.-BlurA(uv, 0)*1.).g*length(GradientA(uv+GradientA(uv, pixelSize*2., vec4(0.,128.,0.,0.), 1)*pixelSize, pixelSize*1.5, vec4(0.,0.,16.,0.), 0)));

    finalColor *= 0.5+(sin(_uvc.x*5.0+TIME*0.1)*_uv.y+cos(_uvc.x*4.7-TIME*0.47)+cos(_uvc.x*11.0+TIME*0.89)*(1.0-_uv.y))/6.0;


    fragColor = vec4(finalColor, 1.0);

    //    fragColor = BlurA(uv, 0); // simple bypass
    //    fragColor = BlurB(uv, 0); // simple bypass
       // fragColor = texture(buffD, _uv); // raw Gaussian pyramid

}



vec4 renderMain(void)
{
  if (PASSINDEX == 0.0){
    //buffA
    vec4 fragColor;
    mainImage0(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 1.0){
    //buffC
    vec4 fragColor;
    mainImage2(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 2.0){
    //Image
    vec4 fragColor;
    mainImage3(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 3.0){
    //Image
    vec4 fragColor;
    mainImage4(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
}
