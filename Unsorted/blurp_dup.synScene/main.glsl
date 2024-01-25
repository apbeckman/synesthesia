vec2 pixelSize = 1./ RENDERSIZE.xy;
vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
float smin(float a, float b, float k) {
    float h = max(k - abs(a-b), 0.) / k;
    return min(a, b) - h*h*h*k*1./6.;
}

float iGlobalTime = TIME*1.0;
int iFrame = int(TIME*60.0);

// vec2 lightPos = vec2(resolution.x*(0.5+0.5*sin(vuTime*2*PI*0.5)),resolution.y*(0.5+0.5*cos(vuTime*2*PI*0.67)));

#define pi2_inv 0.159154943091895335768883763372
#define pi 3.141592653589793238462643383279
		
vec2 complex_mul(vec2 factorA, vec2 factorB)
{
  return vec2( factorA.x*factorB.x - factorA.y*factorB.y, factorA.x*factorB.y + factorA.y*factorB.x);
}

vec2 complex_div(vec2 numerator, vec2 denominator){
  return vec2( numerator.x*denominator.x + numerator.y*denominator.y,
    numerator.y*denominator.x - numerator.x*denominator.y)/
    vec2(denominator.x*denominator.x + denominator.y*denominator.y);
}

float sigmoid(float x)
{
  return 2./(1. + exp2(-x)) - 1.;
}

float border(vec2 domain, float thickness){
  vec2 uv = fract(domain-vec2(0.5));
  uv = min(uv,1.-uv)*2.;
  return clamp(max(uv.x,uv.y)-1.+thickness,0.,1.)/(thickness);
}

float square_mask(vec2 domain){
  return (domain.x <= 1. && domain.x >= 0. && domain.y <= 1. && domain.y >= 0.) ? 1. : 0.;
}

// see https://stackoverflow.com/a/26070411
float atan2(float y, float x)
{
  float s = float(abs(x) > abs(y));
  return mix(pi/2.0 - atan(x,y), atan(y,x), s);
}

vec2 spiralzoom(vec2 domain, vec2 center, float n, float spiral_factor, float zoom_factor, vec2 pos){
  vec2 uv = domain - center;
  float d = length(uv);
  return vec2( atan2(uv.y, uv.x)*n*pi2_inv + d*spiral_factor, -log(d)*zoom_factor) + pos;
}

float mirror(float v){
  return 1.-abs(fract(abs(v*0.5))*2.-1.);
}

vec2 mirror(vec2 v){
  return vec2( mirror(v.x), mirror(v.y) );
}

// HSL to RGB converter code from http://www.gamedev.net/topic/465948-hsl-shader-glsl-code/
float Hue_2_RGB(float v1, float v2, float vH )
{
  float ret;
  if ( vH < 0.0 )
    vH += 1.0;
  if ( vH > 1.0 )
    vH -= 1.0;
  if ( ( 6.0 * vH ) < 1.10 )
    ret = ( v1 + ( v2 - v1 ) * 6.0 * vH );
  else if ( ( 2.0 * vH ) < 1.0 )
    ret = ( v2 );
  else if ( ( 3.0 * vH ) < 2.0 )
    ret = ( v1 + ( v2 - v1 ) * ( ( 2.0 / 3.0 ) - vH ) * 6.0 );
  else
    ret = v1;
  return ret;
}

vec3 hsl2rgb(float H, float S, float L){
  float var_2, var_1, R, G, B;
  if (S == 0.0)
  {
    R = L;
    G = L;
    B = L;
  }
  else
  {
    if ( L < 0.5 )
    {
      var_2 = L * ( 1.0 + S );
    }
    else
    {
      var_2 = ( L + S ) - ( S * L );
    }

    var_1 = 2.0 * L - var_2;

    R = Hue_2_RGB( var_1, var_2, H + ( 1.0 / 3.0 ) );
    G = Hue_2_RGB( var_1, var_2, H );
    B = Hue_2_RGB( var_1, var_2, H - ( 1.0 / 3.0 ) );
  }
  return vec3(R,G,B);
}
		
float lum(vec4 rgb)
{
  return dot(rgb, vec4(0.3, 0.59, 0.11, 0.));
}

vec2 gradient(sampler2D sampler, vec2 uv, vec2 d, vec4 selector)
{
  vec4 dX = texture(sampler, uv + vec2(1.,0.)*d) - texture(sampler, uv - vec2(1.,0.)*d);
  vec4 dY = texture(sampler, uv + vec2(0.,1.)*d) - texture(sampler, uv - vec2(0.,1.)*d);
  return vec2( dot(dX, selector), dot(dY, selector) );
}

vec2 rot90(vec2 vector){
  return vector.yx * vec2(1,-1);
}

float circle_distance(vec2 uv, vec2 pos, float size, float min)
{
  return max( min, 1. - length((uv - pos) * aspect / size) );
}

float smooth_circle(vec2 uv, vec2 pos, float size, float ramp)
{
  return 0.5 + sigmoid( circle_distance(uv, pos, size, -16.) * ramp) * 0.5;
}

vec2 vortex_warp(vec2 uv, vec2 pos, float size, float ramp, vec2 rot)
{
  vec2 pos_correct = 0.5 + (pos - 0.5);
  vec2 rot_uv = pos_correct + complex_mul((uv - pos_correct)*aspect, rot)/aspect;
  float smooth_circle = smooth_circle(uv, pos_correct, size, ramp);
  return mix(uv, rot_uv, smooth_circle);
}

vec2 vortex_pair_warp(vec2 uv, vec2 pos, vec2 vel)
{
  vec2 aspect = vec2(1.,RENDERSIZE.y/RENDERSIZE.x);
  float ramp = 5.;
  float d = 0.2;

  float l = length(vel);
  vec2 p1 = pos;
  vec2 p2 = pos;

  if(l > 0.)
  {
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

vec4 noise()
{
  vec2 uv = gl_FragCoord.xy * pixelSize * RENDERSIZE / vec2(256);
  return texture(colornoise, uv + rnd.xy);
}

float lum(vec3 colIn){
  return dot(colIn, vec3(1.0))/3.0;
}

vec4 advance()
{
  vec2 uv = gl_FragCoord.xy * pixelSize;
  vec4 prev = texture(sampler_prev, uv);
  vec4 blur = texture(sampler_blur, uv);
  vec4 blur2 = texture(sampler_blur2, uv);
  vec4 blur3 = texture(sampler_blur3, uv);
  vec4 noise = noise();

  vec2 brushPosition = vec2((sin(TIME*0.3+syn_BassTime*2*PI*0.05)),(cos(TIME*0.25+syn_BassTime*2*PI*0.04)));
  float stir = manual_stir;
  // kindly borrowed from glassier scene
  if (_mouse.z > 0. && stir_or_paint != 1.0){
    brushPosition = _muvc;
    stir = _mouse.z*0.000125;
  }
  else{
    stir = manual_stir;
  }

  float twitcher = ( pow(syn_BassLevel, 1.3) * auto_stir * 1.25 + stir * 4. ) * stir_strength;

  uv = vortex_pair_warp(uv, brushPosition*0.5+0.5, normalize(brushPosition)*(twitcher)*0.025);

  float scaleModder = mix(1.0, pow(syn_BassLevel*0.46 + syn_MidLevel * 0.46,0.9)*2.0, auto_zoom)*(zoom_in_out*0.6+0.1*zoom_in_out*abs(zoom_in_out));

  vec2 centerModded = _uvc - (center-0.5);

  // translate and scale with a zoom factor controlled by difference of low and high freq volume levels
  uv = center + (uv - center) * (1.0 - vec2(scaleModder*0.25 + shock_forward * 2.) * aspect * pixelSize * 64.);
  uv += centerModded*length(centerModded)*0.05*(fisheye_scale*fisheye_scale)*sign(fisheye_scale);
  uv += centerModded*(0.5+0.5*length(centerModded))*0.05*shock_back;
  uv += translate * pixelSize;
  uv -= vec2(center.yx - uv.yx) * aspect.yx * shear * pixelSize * 256.;

  // uv += _rotate(vec2(1.0,0.0), 3*PI*_fbm(vec3(_uvc*0.5, syn_BassTime*0.5)))*0.01*melter*syn_BassLevel;

  uv += _rotate(vec2(1.0,0.0), 3*PI*_fbm(vec3(_uvc*0.25, syn_BassTime*0.3)))*0.0035*warper*pow(syn_BassLevel, 1.5);
  uv += prev.r*_rotate(vec2(1.0,0.0), 8*PI*pow(_fbm(vec3(_uvc*8.0, 0.3*syn_MidHighTime)), 2.0))*0.005*fritzer*pow(syn_HighHits,1.5);


  prev = texture(sampler_prev, uv); // sample with warped lookup vector
  prev.g = texture(sampler_prev, uv+0.005*_rotate(vec2(0.0,1.0), blur.g*3.33*PI)).g; // sample with warped lookup vector
  float r = prev.r + (blur2.r - blur3.r) * 32. / 256.;  // reaction-diffusion
  float g = prev.g + (blur2.g - blur3.g) * 64. / 256.;  // reaction-diffusion

  if (_mouse.z > 0. && stir_or_paint == 1.0 && (length(_xy - _mouse.xy)) < 140*paint_size){
       r += lum(vec3(smin(PI*2*_noise(r), PI*2*_fbm(r), 0.5)));
       g += lum(vec3(smin(PI*8*_noise(g), PI*8*_fbm(g), 0.5)));
  }
  else{
    r = r;
  }
  r += (noise.x - 0.5) * 32. / 1920.;
  g += (noise.x - 0.5) * 32. / 1920.;

  r = clamp(r, 0., 1.);
  r += media_burn*sin(lum(_loadUserImage())*2*PI);

  r += center_seed*(-1.0+fract(length(_uvc)*10.0)*2.0);

  g = clamp(g, 0., 1.);
  g += media_burn*sin(lum(_loadUserImage())*8*PI);
  g += center_seed*(-1.0+fract(length(_uvc)*10.0)*3.0);


  if(FRAMECOUNT <= 10)
  {
    return noise;
  }
  else
  {
    return vec4(r, g, 0, r);
  }
}


vec3 _grad4(vec3 col1, vec3 col2, vec3 col3, vec3 col4, float mixVal){
    mixVal *= 3.0;
    float mix1 = clamp(mixVal,0.0,1.0);
    float mix2 = clamp(mixVal-1.0, 0.0, 1.0);
    float mix3 = clamp(mixVal-2.0, 0.0, 1.0);
    vec3 firstTwo = mix(mix(col1, col2, mix1), mix(col2, col3, mix2), step(1.0, mixVal));
    return mix(firstTwo, mix(col3, col4, mix3), step(2.0, mixVal));
}

vec4 composite()
{
  vec2 uv = gl_FragCoord.xy / RENDERSIZE.xy;

  vec4 noise = noise();
  vec4 rawPrev = texture(sampler_prev, uv);
  vec4 rawBlur2 = texture(sampler_blur2, uv);
  vec4 rawBlur3 = texture(sampler_blur3, uv);

  float interestingMixer = _grad4(vec3(0.0), 0.5+0.5*vec3(rawPrev.g), 0.5+0.5*vec3(rawPrev.r), vec3(1.0), roil).r;

  rawPrev = mix(rawPrev.rgba, rawPrev.grba, interestingMixer);
  rawBlur2 = mix(rawBlur2.rgba, rawBlur2.grba, interestingMixer);
  rawBlur3 = mix(rawBlur3.rgba, rawBlur3.grba, interestingMixer);

  vec4 colOut = rawBlur3.r*vec4(0);

  vec2 grad1 = gradient(sampler_blur, uv, pixelSize * 4., vec4(1,0,0,0) );
  vec2 grad3 = gradient(sampler_blur3, uv, pixelSize * 4., vec4(1,0,0,0) );
  vec2 grad1B = gradient(sampler_blur, uv, pixelSize * 4., vec4(0,1,0,0) );
  vec2 grad3B = gradient(sampler_blur3, uv, pixelSize * 4., vec4(0,1,1,0) );
  grad1 = mix(grad1, grad1B, interestingMixer);
  grad3 = mix(grad3, grad3B, interestingMixer);

  //vec4 col1 = mix(vec4(0.1,0.0,0.7,0), vec4(1,0.7,0.3,0), color_mode*rawPrev.r);
  vec4 col1 = mix(vec4(0.6,0.5,0.7,0), vec4(1,0.6,0.3,0), color_mode*rawPrev.r);
  //vec4 col2 = mix(vec4(0,1,1,0), vec4(0.3,0,0.4,0), color_mode*rawPrev.r);
  vec4 col2 = mix(vec4(0,1,1,0), vec4(0.3,0,0.4,0), color_mode*rawPrev.r);

  vec2 light_angle_mod = light_angle;
  if (auto_circling>0.5){
    light_angle_mod = vec2(0.5+0.5*sin(smoothTimeC*0.3), 0.5+0.5*cos(.1*smoothTime));
  }

  colOut = mix(colOut, mix(col1, col2, pow(length(uv - light_angle_mod + grad1*1.)*1.5,2.)), mirror(length(grad3)*4.));
  //vec3 flashGradient = _grad4(vec3(0.0,0.0,0.0), vec3(0.9,0.7,0.0), vec3(0.1,0.7,0.9), vec3(1.0), pow(mirror(length(grad3)*4.),2.0)); 
  vec3 flashGradient = _grad4(vec3(0.0,0.0,0.0), vec3(0.9,0.7,0.0), vec3(0.1,0.7,0.9), vec3(1.0), pow(mirror(length(grad3)*4.),2.0)); 

  colOut += vec4(flashGradient, 0.0)*_pulse(length(_uvc)+grad1.r, fract(syn_MidHighTime*0.2)*2.0-0.5, 0.5); 
  // colOut = mix(colOut, vec4(0.0,0.6,0.7,0.0), color_mode*rawPrev.r);

  // float bottoms = pow(_pulse(0.2-grad1.r, 0.0, 0.3)-_pulse(0.2-grad3.r, 0.0, 0.3), 2.0);
  float bottoms = pow(_pulse(0.2-grad1.r, 0.0, 0.3), 2.0);
  vec4 lightningCol = vec4(_grad4(vec3(0.0), vec3(0.0,0.0,0.9), vec3(0.0,0.6,1.0), vec3(0.5,1.0,1.0), bottoms), 0.0);

  colOut += syn_HighHits*lightningCol*_pulse(grad3.r*3.0, fract(TIME*0.3)*1.0-0.25, 0.2);

  // float veryBottom = pow(_pulse(rawBlur2.r-rawBlur3.r, 0.0, 0.1),1);
  // float pulsePos = ((_uv.x+_uv.y)*0.5+veryBottom*0.5)/1.5;
  // float thickness = 0.1;
  // float network = _pulse(pulsePos, fract(TIME*0.2)-thickness*0.5, thickness);
  // colOut += network;
  // vec4 pixCol = vec4(_pixelate(colOut.rgb, 5.0),_pixelate(colOut.a, 5.0));
  // colOut = mix(colOut, pixCol, posterize);

  /*  
  float c1 = circle_distance(uv, vec2(0.25, 0.25 + 0.5 * syn_BassLevel), pixelSize.y * 25., 0.);
  float c2 = circle_distance(uv, vec2(0.75, 0.25 + 0.5 * syn_HighLevel), pixelSize.y * 25., 0.);
  output = mix(output, vec4(1,0,0,0), c1);
  output = mix(output, vec4(0,1,1,0), c2);
  */

  return colOut;
}

vec4 blurH(sampler2D sampler_prev)
{
  vec2 uv = gl_FragCoord.xy * pixelSize;
  float h = pixelSize.x;
  vec4 sum = vec4(0.0);
  // Gaussian factors from https://github.com/mattdesl/lwjgl-basics/wiki/ShaderLesson5
  // horizontal blur fragment shader
  sum += texture(sampler_prev, vec2(-4,0) * pixelSize + uv ) * 0.0162162162;
  sum += texture(sampler_prev, vec2(-3,0) * pixelSize + uv ) * 0.0540540541;
  sum += texture(sampler_prev, vec2(-2,0) * pixelSize + uv ) * 0.1216216216;
  sum += texture(sampler_prev, vec2(-1,0) * pixelSize + uv ) * 0.1945945946;
  sum += texture(sampler_prev, vec2( 0,0) * pixelSize + uv ) * 0.2270270270;
  sum += texture(sampler_prev, vec2( 1,0) * pixelSize + uv ) * 0.1945945946;
  sum += texture(sampler_prev, vec2( 2,0) * pixelSize + uv ) * 0.1216216216;
  sum += texture(sampler_prev, vec2( 3,0) * pixelSize + uv ) * 0.0540540541;
  sum += texture(sampler_prev, vec2( 4,0) * pixelSize + uv ) * 0.0162162162;
  return sum;
}

vec4 blurV(sampler2D sampler_blurH)
{
  vec2 uv = gl_FragCoord.xy * pixelSize;
  float v = pixelSize.y;
  vec4 sum = vec4(0.0);
  sum += texture(sampler_blurH, vec2(0,-4) * pixelSize + uv ) * 0.0162162162;
  sum += texture(sampler_blurH, vec2(0,-3) * pixelSize + uv ) * 0.0540540541;
  sum += texture(sampler_blurH, vec2(0,-2) * pixelSize + uv ) * 0.1216216216;
  sum += texture(sampler_blurH, vec2(0,-1) * pixelSize + uv ) * 0.1945945946;
  sum += texture(sampler_blurH, vec2(0, 0) * pixelSize + uv ) * 0.2270270270;
  sum += texture(sampler_blurH, vec2(0, 1) * pixelSize + uv ) * 0.1945945946;
  sum += texture(sampler_blurH, vec2(0, 2) * pixelSize + uv ) * 0.1216216216;
  sum += texture(sampler_blurH, vec2(0, 3) * pixelSize + uv ) * 0.0540540541;
  sum += texture(sampler_blurH, vec2(0, 4) * pixelSize + uv ) * 0.0162162162;
  return sum;
}

vec4 renderMain(void)
{
  if      (PASSINDEX == 0.0)  return advance();
  else if (PASSINDEX == 1.0)  return texture(sampler_prev_B, _uv);
  else if (PASSINDEX == 2.0)  return blurH(sampler_prev);
  else if (PASSINDEX == 3.0)  return blurV(sampler_blurH);
  else if (PASSINDEX == 4.0)  return blurH(sampler_blur);
  else if (PASSINDEX == 5.0)  return blurV(sampler_blur2H);
  else if (PASSINDEX == 6.0)  return blurH(sampler_blur2);
  else if (PASSINDEX == 7.0)  return blurV(sampler_blur3H);
  else if (PASSINDEX == 8.0) return composite();
}
