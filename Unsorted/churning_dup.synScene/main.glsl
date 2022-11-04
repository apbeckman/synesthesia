float vuTimeUncorrected = syn_Time;
float vuFadeToBlack = syn_FadeInOut;
float highHits = syn_HighHits;
float highAccum = syn_HighLevel;
float bassAccum = syn_BassLevel;
float bassHits = syn_BassHits;
float switchOnBeat = syn_ToggleOnBeat;
float randomizerBeat = syn_RandomOnBeat;

float vuTime = vuTimeUncorrected*60000;
float time = vuTime;
float beatTime = syn_BeatTime*60000;

vec2 resolution = RENDERSIZE;

// const float PI = 3.141592;

vec4 sigmoid(vec4 x){
  return x/(1.+abs(x));
}

void firstPass( out vec4 fragColor, in vec2 fragCoord )
{
  // float t0 = cos(.1*time);
  // float t1 = cos(.17*time);
  // float t2 = cos(.117*time);
  // float t3 = cos(.1717*time);
  // float t4 = cos(.11717*time);

  vec2 dirOfPush = normalize(vec2(sin(direction*2*PI+syn_BPMTri2), cos(direction*2*PI+randomizerBeat)));

  if (auto_push<0.5){
    dirOfPush = vec2(1.0,1.0)*manual_push;
  }

  vec2 tAll = (syn_Level*syn_Level)*2.0*dirOfPush*(0.8+syn_HighHits*0.8);

  if (lift > 0.05){
    dirOfPush = vec2(-1.0,0.0);
    tAll = lift*10.0*dirOfPush;
  }

  if (swirl > 0.05){
    dirOfPush = normalize(_uvc);
    tAll = swirl*40.0*dirOfPush;
  }


  vec2 t0 = _rotate(tAll, time);
  vec2 t1 = _rotate(tAll, -time);
  vec2 t2 = _rotate(tAll, PI/2);
  vec2 t3 = _rotate(tAll, -PI/2);
  vec2 t4 = -tAll;

  vec2 wp = 4.0/resolution.xy;
  float wg =  -5.0;

  vec2 uv = fragCoord.xy / resolution.xy;
  uv += _rotate(vec2(1.0,0.0),length(_uvc))*pow(dive_in,2.0);

  vec4 c = -1.+2.*texture(firstFB, uv);
  c = -1.+2.*texture(firstFB, uv + t0*wp*c.rg);
  c = -1.+2.*texture(firstFB, uv - t1*wp*c.ba);
  c = -1.+2.*texture(firstFB, uv + t2*wp*c.ar);
  c = -1.+2.*texture(firstFB, uv - t3*wp*c.gb);

  vec4 lp = -1.+2.*texture(secondFB, uv);

  // lp *= equalizer3;

  // uv = uv*(1.0+fract(length(_uvc)*equalizer3));

  vec4 rates = vec4(1.,.11,.0111,.001111);

  vec4 phi = rates*time;
  phi = vec4(highHits);

  // vec4 g = vec4(dot(uv,vec2(1.,0.)+wg*c.rb), dot(uv, vec2(.3,0.8)+wg*c.ga), dot(uv,vec2(.7,.7)+wg* c.bg), dot(uv, vec2(0.,1.)+wg*c.ar));
  // g = sin(phi+2.*PI*g + wg*c.gbar);
  vec4 g = vec4(
                dot(fract(uv-c.rg*broad_strokes), vec2(1.,1.)),
                dot(fract(uv+c.ba*broad_strokes), vec2(-2.,0.)),
                dot(fract(uv+c.gr*broad_strokes), vec2(0.,2.)),
                dot(fract(uv+c.ab*broad_strokes), vec2(1.,-1.))
                );
  g = sin(2.0*PI*g + wg*c.gbar + phi);

  g = sigmoid(2.*g.rgba-g.gbar);
  g *= vec4(1., .57, .27, .77)*(1.0+syn_HighHits);


  vec4 c_new = .5*c + g + .25*sigmoid(1.*(c-lp));
  c_new = c_new.gbar;

  vec4 c_diff = c_new - c;
  c_diff = .04*mix(c_diff, c_diff.gbra - c_diff.brga, .3*length(t4)+.4);

  c_diff += (plasma+.05)*normalize(c_new-dot(c_new, vec4(.25)));

  float presence;
  if (syn_HighPresence>syn_BassPresence){
    presence = syn_HighPresence;
  } else if (syn_HighPresence<syn_BassPresence){
    presence = -syn_BassPresence;
  }

  // if (s2 > 0.5){
  //   c_new = (1.0-s3)+(s3)*clamp(c+c_diff,-1.0, 1.0);
  // } else {
    c_new = .5+.5*clamp(c+c_diff*(clamp(presence,-1.0,1.0)+0.05),-1.,1.);
  // }

  // c_new = .5+.5*clamp(c+c_diff*(clamp(presence,-1.0,1.0)+0.05),-1.,1.);

  if (turbulence > 0.5){
    c_new = c_new - c_diff;
  }

  if (_exists(syn_UserImage)){
    vec4 image = _loadUserImage();
    c_new = mix(image,c_new,1.0-img_mix);
  }

  // c_new *= vec4(equalizer1, equalizer2, equalizer3, 1.0);

  fragColor = c_new;
}


//Second
void secondPass( out vec4 fragColor, in vec2 fragCoord )
{
  vec2 uv = fragCoord.xy / resolution.xy;
  vec2 d = 1./resolution.xy;
  vec4 c = .2*(texture(firstFB, uv)
               + texture(secondFB, uv+d*vec2(1.,0.))
               + texture(secondFB, uv+d*vec2(-1.,0.))
               + texture(secondFB, uv+d*vec2(0.,-1.))
               + texture(secondFB, uv+d*vec2(0.,1.))
               );
  fragColor = c;
}




//Third
vec3 sigmoid(vec3 x){
  return x/(1.+abs(x));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
  vec2 uv = fragCoord.xy / resolution.xy;
  vec4 s = texture(firstFB, uv)*2.-1.;
  vec4 s2 = texture(secondFB, uv)*2.-1.;
  s = mix(s2, s, highHits);
  vec3 c = .5+.5*sigmoid(2.*(s.rgb - .25*s.bar));
  c = c + .05 + .15*normalize(c-dot(c,vec3(.3333)));
  // fragColor = vec4(c.rgb*vec3(1.0,1.0,0.5+highAccum*_fbm(_uvc-vuTime)),1.);
  // if (colRegSel == 0.0){

  // }

  float mixer = fract(syn_Time*0.2);
  float reg;
  float timer = modf(mod(syn_Time*0.2,3.0), reg);

  vec2 changer = vec2(1.0-mixer, mixer);
  vec3 v = vec3(0.0);
  if (reg == 0.0){
    v = vec3(changer, 0.0);
  }
  if (reg == 1.0){
    v = vec3(0.0, changer);
  }
  if (reg == 2.0){
    v = vec3(changer.y, 0.0, changer.x);
  }

  if (color_phasing > 0.5){
    c = vec3(c.r*v.r + c.g*v.g + c.b*v.b, c.g*v.r + c.b*v.g + c.r*v.b, c.b*v.r + c.r*v.g + c.g*v.b);
  }

  if (colRegSel == 0.1){
    c = c.bgr-c.rrr+vec3(0.0,0.0,c.g*0.2);
    c *= 1.5;
    c *= vec3(1.0, 1.0, 1.0+syn_HighHits*0.5);
  }
  else if (colRegSel == 0.2){
    c = (c.gbr+c.brb)*0.5;
    c= pow(c,vec3(3.0));
    c *= vec3(1.0, 1.0+syn_HighHits*0.5, 1.0+syn_HighHits)*1.2;
  }

  fragColor = vec4(c,1.);
}

float SCurve (float x) {
        x = x * 2.0 - 1.0;
        return -x * abs(x) * 0.5 + x + 0.5;

        // return dot(vec3(-x, 2.0, 1.0 ),vec3(abs(x), x, 1.0)) * 0.5; // possibly faster version

}

vec4 BlurV (sampler2D source, vec2 size, vec2 uv, float radius) {

    if (radius >= 1.0)
    {
        vec4 A = vec4(0.0);
        vec4 C = vec4(0.0);

        float height = 1.0 / size.y;

        float divisor = 0.0;
        float weight = 0.0;

        float radiusMultiplier = 1.0 / radius;

        for (float y = -radius; y <= radius; y++)
        {
            A = texture(source, uv + vec2(0.0, y * height));

                weight = SCurve(1.0 - (abs(y) * radiusMultiplier));

                C += A * weight;

            divisor += weight;
        }

        return vec4(C.r / divisor, C.g / divisor, C.b / divisor, 1.0);
    }

    return texture(source, uv);
}

vec4 BlurH (sampler2D source, vec2 size, vec2 uv, float radius) {

    if (radius >= 1.0)
    {
        vec4 A = vec4(0.0);
        vec4 C = vec4(0.0);

        float width = 1.0 / size.x;

        float divisor = 0.0;
        float weight = 0.0;

        float radiusMultiplier = 1.0 / radius;

        // Hardcoded for radius 20 (normally we input the radius
        // in there), needs to be literal here

        for (float x = -radius; x <= radius; x++)
        {
            A = texture(source, uv + vec2(x * width, 0.0));

                weight = SCurve(1.0 - (abs(x) * radiusMultiplier));

                C += A * weight;

            divisor += weight;
        }

        return vec4(C.r / divisor, C.g / divisor, C.b / divisor, 1.0);
    }

    return texture(source, uv);
}

float hZone = (1 - abs(mod(syn_BassTime*0.05, (2)) - 1))-0.5;
float tiltShiftZone = (1.0-smoothstep(0.0+hZone, 0.5+hZone, _uv.y))+(smoothstep(0.5+hZone, 1.0+hZone, _uv.y));

vec4 vertBlurPass(){
    return BlurV(vertBlur, RENDERSIZE, _uv, tiltShiftZone*20.0*blur_on_fps);
}

vec4 horBlurPass(){
    return BlurH(horBlur, RENDERSIZE, _uv, tiltShiftZone*20.0*blur_on_fps);
}


vec4 renderMain () {
  vec4 fragColor = vec4(0.0);
  float blurTime = (1.0+0.5*sin(syn_Time*0.1)+0.5*sin(syn_Time*0.0177))*0.5;
  if (PASSINDEX == 0.0){
    firstPass(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 1.0){
    secondPass(fragColor, gl_FragCoord.xy);
    return fragColor;
  }
  else if (PASSINDEX == 2.0){
    mainImage(fragColor, gl_FragCoord.xy);
    return fragColor;
  }   
  else if (PASSINDEX == 3.0){
    return vertBlurPass();
  }
  else if (PASSINDEX == 4.0){
    return horBlurPass()-tiltShiftZone*0.05;
  }

  return vec4(1.0, 0.0, 0.0, 1.0);
}
