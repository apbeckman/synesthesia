float bassHits = syn_BassHits;
float highHits = syn_HighHits;
float vuTimeUncorrected = syn_Time;
float bassTimeUncorrected = syn_BassTime;
float randomizerBeat = syn_RandomOnBeat;
float midsAccum = syn_MidLevel;
float beatTimeUncorrected = syn_BeatTime;
float onBeat = syn_OnBeat;
float switchOnBeat = syn_ToggleOnBeat;
float highsAccum = syn_HighLevel;
float bassAccum = syn_BassLevel;

vec3 iResolution = vec3(RENDERSIZE.x, RENDERSIZE.y, 0.0);
float iGlobalTime = TIME*1.0;
int iFrame = int(TIME*60.0);
float time = TIME;
vec2 resolution = RENDERSIZE;
float vuTime = vuTimeUncorrected*6*0.01 + time * 0.05;
float bassTime = syn_BassTime*0.09;
float beatTime = beatTimeUncorrected*6;

vec2 lightPos = vec2(resolution.x*(0.5+0.5*sin(bassTime*2*PI*0.18)),resolution.y*(0.5+0.5*cos(bassTime*2*PI*0.224)));

vec4 iMouse = vec4(lightPos,bassAccum*2.0,1.1);

//General Functions

#define pi2_inv 0.159154943091895335768883763372

float sdRoundBox(vec2 p, vec2 b, float r)
{
  vec2 q = abs(p) - b;
  vec2 m = vec2( min(q.x,q.y), max(q.x,q.y) );
  float d = (m.x > 0.0) ? length(q) : m.y;
  return d - r;
}

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

vec2 mouseDelta(){
  vec2 pixelSize = 1. / iResolution.xy;
  float eighth = 1./8.;
  vec4 oldMouse = texture(buffD, vec2(7.5 * eighth, 2.5 * eighth));
  vec4 nowMouse = vec4(iMouse.xy / iResolution.xy, iMouse.zw / iResolution.xy);
  if(oldMouse.z > pixelSize.x && oldMouse.w > pixelSize.y &&
     nowMouse.z > pixelSize.x && nowMouse.w > pixelSize.y)
  {
    return nowMouse.xy - oldMouse.xy;
  }
  return vec2(0.);
}

float border(vec2 domain, float thickness){
  vec2 uv = fract(domain-vec2(0.5));
  uv = min(uv,1.-uv)*2.;
  return clamp(max(uv.x,uv.y)-1.+thickness,0.,1.)/(thickness);
}

vec2 rot90(vec2 vector){
  return vector.yx*vec2(1,-1);
}

vec2 wrap_flip(vec2 uv){
  return vec2(1.)-abs(fract(uv*(.5-simplify*0.4))*2.-1.);
}

vec2 complex_mul(vec2 factorA, vec2 factorB){
  return vec2( factorA.x*factorB.x - factorA.y*factorB.y, factorA.x*factorB.y + factorA.y*factorB.x);
}

vec2 rotozoom(vec2 uv, float ang, float zoom, vec2 aspect){
  vec2 rot = vec2(cos(ang), sin(ang))*zoom;
  return 0.5 + complex_mul((uv - 0.5)*aspect, rot)/aspect;
}

vec2 spiralzoom(vec2 domain, vec2 center, float n, float spiral_factor, float zoom_factor, vec2 pos){
  vec2 uv = domain - center;
  float d = length(uv);
  return vec2( atan(uv.y, uv.x)*n*pi2_inv + d*spiral_factor, -log(d)*zoom_factor) + pos;
}

vec2 complex_div(vec2 numerator, vec2 denominator){
  return vec2( numerator.x*denominator.x + numerator.y*denominator.y,
              numerator.y*denominator.x - numerator.x*denominator.y)/
  vec2(denominator.x*denominator.x + denominator.y*denominator.y);
}

float circle(vec2 uv, vec2 aspect, float scale){
  return clamp( 1. - length((uv-0.5)*aspect*scale), 0., 1.);
}

float sigmoid(float x) {
  return 2./(1. + exp2(-x)) - 1.;
}

float smoothcircle(vec2 uv, vec2 aspect, float radius, float ramp){
  return 0.5 - sigmoid( ( length( (uv - 0.5) * aspect) - radius) * ramp) * 0.5;
}

float conetip(vec2 uv, vec2 pos, float size, float min)
{
  vec2 aspect = vec2(1.,iResolution.y/iResolution.x);
  return max( min, 1. - length((uv - pos) * aspect / size) );
}

float warpf1lter(vec2 uv, vec2 pos, float size, float ramp)
{
  return 0.5 + sigmoid( conetip(uv, pos, size, -16.) * ramp) * 0.5;
}

vec2 vortex_warp(vec2 uv, vec2 pos, float size, float ramp, vec2 rot)
{
  vec2 aspect = vec2(1.,iResolution.y/iResolution.x);

  vec2 pos_correct = 0.5 + (pos - 0.5);
  vec2 rot_uv = pos_correct + complex_mul((uv - pos_correct)*aspect, rot)/aspect;
  float f1lter = warpf1lter(uv, pos_correct, size, ramp);
  return mix(uv, rot_uv, f1lter);
}

vec2 vortex_pair_warp(vec2 uv, vec2 pos, vec2 vel)
{
  vec2 aspect = vec2(1.,iResolution.y/iResolution.x);
  float ramp = 4.;

  float d = 0.125;

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

  vec4 loadUserImageZoomed(vec2 offset, float zoom) {
    vec2 uv = _correctImageCoords(textureSize(syn_UserImage, 0));
    uv *= zoom;
    uv -= 1/zoom;
    
    vec4 tempCol = _invertImage(_imageToMask(texture(syn_UserImage, uv+offset)));

    if (uv.x < 0.0 || uv.y < 0.0 || uv.y > 1.0 || uv.x > 1.0) {
      return _invertImage(vec4(0.0));
    }
    
    return tempCol;
  }


// ****************** PASS 0 ***********************
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
  return vec2( dot(dX, selector), dot(dY, selector) );
}

void mainImage0(out vec4 fragColor, in vec2 fragCoord)
{
  vec2 uv = fragCoord.xy / iResolution.xy;

  vec2 aspect = vec2(1.,iResolution.y/iResolution.x);
  vec2 pixelSize = 1. / iResolution.xy;

  float center_f1lter = smoothcircle(uv, aspect, 0.1, 32.*(1.4-syn_HighPresence*0.4-0.4*syn_MidPresence));

  float ang = (bassTime*0.06+syn_BPMTwitcher*0.06*beat_rotate)/pi2_inv;
  float scale = mix(3., 18., center_f1lter);

  vec2 rot = vec2(cos(ang), sin(ang));
  vec2 offset = vec2(pow(sandstorm,3.0)*1000);

  vec2 uv_carpet = 0.5 + complex_mul((uv - 0.5)*aspect, rot*scale)/aspect;
  uv_carpet = wrap_flip(uv_carpet + vec2(1.0)-_uvc*offset*0.1);


  float carpet_circle = circle(uv_carpet, aspect, 4.);


  // blue is fed by a carpet-mirrored circle and a carpet-transformed texture feedback of itself
  fragColor.b = BlurA(uv_carpet, 0).b * mix(0.75, 1., center_f1lter) + mix(0., 0.5, center_f1lter*highHits);
  fragColor.b = mix(fragColor.z, 1., carpet_circle);

  // green is pretty much like blue but with a steeper fall-off, and heavier zoom

  scale = mix(4., 64., center_f1lter);
  vec2 uv_galaxies = 0.5 + complex_mul((uv - 0.5)*aspect, rot*scale)/aspect;
  uv_galaxies = wrap_flip(uv_galaxies + offset);

  float bleedAmt = 0.5;
  fragColor.g = pow(texture(buffA, _uv).g*(0.5+highHits*0.8)+BlurA(uv_galaxies, 0).b*bleedAmt,2.0);

  fragColor.g = mix(fragColor.g, 0.2, center_f1lter);

  // red is fed by green
  vec2 uv_zoom_mix = 0.5 + (uv - 0.5)* mix(0.985, 1.0, center_f1lter);

  fragColor.r = BlurA(uv_zoom_mix, 0).x*mix(0.925, 0.0, center_f1lter) - mix(0.02, 0.35, center_f1lter);
  fragColor.r = mix(fragColor.x, 1.,(1.-center_f1lter) * fragColor.y*(1.0-highHits));

  if (_exists(syn_UserImage)){
    vec4 tempCol = loadUserImageZoomed(vec2(0.0),2.0);
    tempCol *= pow(center_f1lter,2.0);
    // fragColor = mix(fragColor+tempCol, fragColor*tempCol, syn_BassLevel);
    fragColor += tempCol*(syn_BassLevel*0.2+0.4*syn_BassPresence+0.3);
  }

  fragColor = clamp(fragColor, 0., 1.);

}


// *************** PASS 1 *****************
//NOT USED, placeholder because shadertoy had it


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

void mainImage2(out vec4 fragColor, in vec2 fragCoord)
{
  vec2 uv = fragCoord.xy / iResolution.xy;
  vec2 uv_half = fract(uv*2.);
  fragColor = blur_horizontal(buffA, uv_half, 1.);
}


// *************** PASS 3 *****************
// vertical blur (second pass)

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


vec2 pattern(vec2 p)
{
  p = fract(p);
  float r = 10.123;
  float v = 0.0, g = 0.0;
  r = fract(r * 9184.928);
  float cp, d;

  d = p.x;
  g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 1000.0);
  d = p.y;
  g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 1000.0);
  d = p.x - 1.0;
  g += pow(clamp(3.0 - abs(d), 0.0, 1.0), 1000.0);
  d = p.y - 1.0;
  g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 10000.0);

  const int iter = 12;
  for(int i = 0; i < iter; i ++)
  {
    cp = 0.5 + (r - 0.5) * 0.9;
    d = p.x - cp;
    g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 200.0);
    if(d > 0.0) {
      r = fract(r * 4829.013);
      p.x = (p.x - cp) / (1.0 - cp);
      v += 1.0;
    }
    else {
      r = fract(r * 1239.528);
      p.x = p.x / cp;
    }
    p = p.yx;
  }
  v /= float(iter);
  return vec2(g, v);
}

// *************** PASS 4 *****************
void mainImage4(out vec4 fragColor, in vec2 fragCoord)
{
  vec2 uv = fragCoord.xy / iResolution.xy;
  vec2 pixelSize = 1. / iResolution.xy;
  vec2 aspect = vec2(1.,iResolution.y/iResolution.x);


  float flowVel = 0.01;
  vec2 flowVector = vec2(sin((onBeat)*PI*2.0+beatTime), cos((onBeat)*PI*2.0+beatTime))*flowVel;

  vec2 grad1_z = GradientA(uv, pixelSize*2., vec4(0,0,0.15,0), 1);
  float center_f1lter = smoothcircle(uv-grad1_z*10.0, aspect, 0.1, 32.);

  vec2 tempVec2 = wrap_flip(uv+grad1_z*highHits);
  vec2 grad0_z = GradientA(tempVec2, pixelSize, vec4(0,0,8,0), 0);


  // float cReg0 = -0.001;
  // float cReg1 = 1.0;
  // float cReg2 = -0.001;
  // float cReg3 = -0.001;
  // float dark = 0.0;
  // float light = 1.0;

  vec3 baseCol = vec3(1.4,1.2,0.2);
  vec3 galaxyCol = vec3(0.0, 0.85, 0.85);
  vec3 subAddColor = vec3(0.50,0.70,0.75);
  vec3 bleed = vec3(0.6,0.8,0.8);
  vec3 fineStructure = vec3(0.282, 0.380, 0.209);
  vec3 depthGalaxies = vec3(0.9);

  // *** Color Regime 1 ***
  baseCol = mix(baseCol, vec3(0.45,0.75,1.0), cReg1);
  galaxyCol = mix(galaxyCol, vec3(1.0,0.7,0.2), cReg1);
  subAddColor = mix(subAddColor, vec3(0.0,0.0,0), cReg1);
  bleed = mix(bleed, vec3(1,1,0), cReg1);
  fineStructure = mix(fineStructure, vec3(0.0,0.0,0.0), cReg1);
  depthGalaxies = mix(depthGalaxies, vec3(0.2,0.2,0.8), cReg1);

  // *** Color Regime 2 ***
  baseCol = mix(baseCol, vec3(0.0), cReg2);
  galaxyCol = mix(galaxyCol, vec3(0.45,0.1,1.0), cReg2);
  subAddColor = mix(subAddColor, vec3(0.25,0.55,1.0), cReg2);
  bleed = mix(bleed, vec3(0.25,0.55,1.0), cReg2);
  fineStructure = mix(fineStructure, vec3(1.0,1.0,1.0), cReg2);
  depthGalaxies = mix(depthGalaxies, vec3(0.0,0.3,0.8), cReg2);

  // *** Color Regime 3 ***
  baseCol = mix(baseCol, vec3(0.0,0.0,0.0), cReg3);
  galaxyCol = mix(galaxyCol, vec3(0.9,0.5,0.0), cReg3);
  subAddColor = mix(subAddColor, vec3(0.0,0.5,0.0), cReg3);
  bleed = mix(bleed, vec3(1.0,0.8,0.0), cReg3);
  fineStructure = mix(fineStructure, vec3(1.0,1.0,1.0), cReg3);
  depthGalaxies = mix(depthGalaxies, vec3(1.0,0.7,0.00), cReg3);

  galaxyCol *= (0.5+0.5*syn_BassPresence);
  // bleed *= center_f1lter;

  vec3 finalColor = vec3(0.0);

  finalColor = mix(finalColor, baseCol, 1.-BlurA(uv + grad0_z*0., 0).x*10.0); //General base background

  finalColor = mix(finalColor, galaxyCol, BlurA(uv - grad1_z*highHits, 0).b); // Galaxies/Neurons

  finalColor = finalColor+mix(vec3(-1), vec3(0.5,0.2,0.5), dark)*subAddColor*length(grad0_z); //Subtract or Add

  finalColor -= dot(grad0_z,vec2(1.0))*syn_HighHits*flashing;

  finalColor = mix(finalColor, bleed, BlurA(uv + grad1_z*5.0, 0).x); // Darkness bleeding away from teal

  vec3 backPat = vec3(1.) * BlurA(uv - GradientA(uv, pixelSize*2., vec4(0,0,2,0), 1), 0).z*(1.-center_f1lter);
  

  vec2 sqPat = pattern(_toPolar((uv-0.5)*1.0*vec2(resolution.x/resolution.y,1.0))+backPat.x*0.05+ flowVector);
  float pulseTime = fract(length(_uvc)-syn_BPMTri2/syn_BPMConfidence + sqPat.x*0.01)*0.5;
  float pulse = smoothstep(0.0,0.2,pulseTime)-smoothstep(0.2,0.3,pulseTime);

  float floatSqPat = step(1.4, sqPat.x)*(1.-center_f1lter)*pulse*radial_grid;
  finalColor = mix(finalColor, fineStructure*floatSqPat,  BlurA(uv, 0).g*floatSqPat - clamp(center_f1lter*1.15-0.25, 0., 1.)); //Fine Structure

  finalColor += floatSqPat*syn_Presence*1.0;

  finalColor = mix(finalColor, depthGalaxies, BlurA(uv - GradientA(uv, pixelSize*2.5, vec4(0,0,3.5,0), 1), 0).z*(1.-center_f1lter)*(1.25-center_f1lter*1.)-0.5); // Depth to Galaxies

  // finalColor = mix(finalColor*0.9, finalColor, center_f1lter*bassAccum);

  float vignetteMod = sdRoundBox(2*(_uv-0.5), vec2(0.6), 0.1);//Magic num 0.2 for nice edges
  finalColor = finalColor*smoothstep(0.5,1.0,1.0-vignetteMod);//Generate Vignette Modifier

  float pulseTime2 = fract(length(_uvc)-syn_BPMTri2/syn_BPMConfidence + pow(length(clamp(finalColor,0.0,1.0)),10.0)*0.01)*0.5;
  float pulse2 = smoothstep(0.0,0.2,pulseTime2)-smoothstep(0.2,0.28,pulseTime2);

  finalColor = mix(finalColor, finalColor*pulse2, pulsate);

  fragColor = vec4(finalColor, 1.0);

}



vec2 getNormPos() {
  return vec2(-1+2*(gl_FragCoord.x/resolution.x),-1+2*(gl_FragCoord.y/resolution.y));
}
vec2 getNormPosBotLeft() {
  return vec2(gl_FragCoord.xy/resolution.xy);
}

vec4 renderMain(void)
{
  if (PASSINDEX == 0.0){
    //buffA
    vec4 fragColor;
    mainImage0(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  // else if (PASSINDEX == 1.0){
  //   //buffB
  //   vec4 fragColor;
  //   mainImage1(fragColor, gl_FragCoord.xy);
  //   return fragColor;
  // }
  else if (PASSINDEX == 1.0){
    //buffC
    vec4 fragColor;
    mainImage2(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 2.0){
    //buffD
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
