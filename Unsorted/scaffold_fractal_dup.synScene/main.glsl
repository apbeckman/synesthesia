float time = TIME;
vec2 mouse = vec2(0.5);
vec2 resolution = RENDERSIZE;

// a raymarching experiment by kabuto
// dmt fork by s1e
// controls, modifications, and audio reactivity by Meebs
// found on glslsandbox.com


const int MAXITER = 25;

vec3 field(vec3 p, float dime) {
    p *= 0.09;
    float f = 0.09;

    for (int i = 0; i < int(iters); i++) {
        p = p.yzx*mat3( .8*scale_xy.x,             .6*scale_xy.y,                             0,
                      -0.6,                        0.8+0.01*sin(dime/120.0),       0, 
                       0.0 + 0.042*sin(dime/25.0), 0.0033 + 0.0012*cos(dime/75.0), 1);
        p += vec3(0.2, cos(TIME/45.0)*0.1 + 0.3, cos(TIME/90.0)*0.1 + 0.3)*float(i);
        p = abs(fract(p)-.5);
        p *= 2.0;
        f *= 2.0;
    }
    p *= p;
    return sqrt(p+p.yzx)/f-.002;
}

float pulse(float variable, float center, float width, float bias){
  return (smoothstep(center-width, center, variable)-(smoothstep(center+bias, center+width+bias, variable)));
}

vec3 rotate_y(vec3 v, float angle)
{
  float ca = cos(angle); float sa = sin(angle);
  return v*mat3(
                +ca, +.0, -sa,
                +.0,+1.0, +.0,
                +sa, +.0, +ca);
}

vec3 rotate_x(vec3 v, float angle)
{
  float ca = cos(angle); float sa = sin(angle);
  return v*mat3(
                +1.0, +.0, +.0,
                +.0, +ca, -sa,
                +.0, +sa, +ca);
}

vec3 opRep( vec3 p, vec3 c )
{
    vec3 q = mod(p,c)-0.5*c;
    return q;
}

vec4 renderMain( void ) {
    float dime = syn_CurvedTime*10.0;
    vec3 dir = normalize(vec3((_xy-resolution*.5)*(perspective+rot_on*4.0)/resolution.x,0.66));
    dir *= vec3(1.0,1.0,z_depth);
    dir = rotate_y(dir, PI*0.5*down_rot_scr);
    dir = rotate_x(dir, PI*0.5*right_rot_scr);
    // dir = mix(diry, dirx, 0.5);

    vec3 pos = vec3(0.0, syn_Presence, dime/16.0);
    vec3 color = vec3(0);
    float hit = 0.0;
    for (int i = 0; i < MAXITER; i++) {
        vec3 f2 = field(pos, dime);
        float f = min(min(f2.x,f2.y),f2.z);
        if (f < 0.0002*100.0){
            hit = 1.0;
        }
        pos += dir*f;
        pos.xy = _rotate(pos.xy, (pos.z-dime/16.0)*pow(rot_on,2.0)*rot_direction);
        color += float(MAXITER-i)/(f2+.000025);
    }
    // vec3 color3 = vec3(1.-1.0/(color*(0.00106)));
    // color3 = color3 * color3 * color3 * color3;

    color = _gamma(vec4(color,1.0), 1.0-pow(rot_on,0.25)*0.25).rgb;

    vec3 posHit = pos*hit;

    float divisor = 0.000106;
    divisor *= 4.5-iters;

    vec3 data = color*divisor;
    vec3 barsData = pow(color*divisor,vec3(10.0));
    vec3 fogData = clamp(data-barsData,0.0,1.0);

    vec3 fogCol = vec3(0.0);
    fogCol += clamp(fogData.r,0.0,1.0)*col3;
    fogCol += fogData.g*col1;
    fogCol += fogData.b*col2;

    // float redPulse = pulse((pos.z-dime/16.0)*2.0, fract(1.0-syn_HighHits)*5.0, syn_HighPresence+0.1, 0.9);
    // float barPulse = smoothstep(syn_BassLevel,0.0,1.0);

    vec3 barCol = vec3(0.0);
    barCol += clamp(barsData.r,0.0,1.0)*col3;
    barCol += clamp(barsData.g,0.0,1.0)*col1;
    barCol += clamp(barsData.b,0.0,1.0)*col2;

    barCol *= bar_on;

    // float barPulse = pulse(posHit.z-dime/16.0,syn_BassPresence+0.2,0.0,0.0);
    // barCol += col3*hit*barPulse;

    // float pulseAmt = _pulse(-0.5+mod(-dime/16.0+pos.z,modFreq)*(1.0/modFreq),fract(syn_BPMTri2/syn_BPMConfidence)*1.0, width);
    float pulseAmt = pulse((posHit.z-dime/16.0), 3.0+fract(smooth_hightime*0.05)*15.0, 0.1, 0.9)*syn_HighHits;

    vec3 finalCol = clamp((fogCol+barCol),0.0,1.0)*(1.0+pulseAmt*10.0);

    // finalCol *= dot(_loadUserImage().rgb,vec3(1.0));

    finalCol.rgb = _gamma(vec4(finalCol, 1.0), 0.5).rgb*(1.5+rot_on*2.0);


    return vec4(finalCol, 1.0);
}