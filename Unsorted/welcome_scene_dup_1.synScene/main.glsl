
vec2 resolution = RENDERSIZE;
float time = syn_Time / 6.0;

vec2 quadrantToPolar(vec2 uv){
  float r = length(uv);
  float theta = 0.0;
  float real = uv.x;
  float imag = uv.y;
  if (real >  0 && imag >= 0) theta = atan(imag/real);
  if (real <= 0 && imag >  0) theta = -1*atan(real/imag)+PI/2;
  if (real <  0 && imag <= 0) theta = atan(imag/real)+PI;
  if (real >= 0 && imag <  0) theta = -1*atan(real/imag)+3*PI/2;
  theta = theta/(2*PI);
  return vec2(r, theta);
}

struct ring
{
  float index;
  float radOut;
  float radIn;
  float squishVal;
  vec2 trans;
};

vec3 colRings[7];

ring ring1 = ring(1.0, 0.445, 0.405, 0.92, vec2(0.0,-0.01));
ring ring2 = ring(2.0, 0.38, 0.34, 0.925, vec2(0.0,-0.01));
ring ring3 = ring(3.0, 0.33, 0.29, 0.95, vec2(0.0,-0.005));
ring ring4 = ring(4.0, 0.281, 0.245, 1.0, vec2(0.002,-0.000));
ring ring5 = ring(5.0, 0.245, 0.203, 1.0, vec2(0.002,-0.000));
ring ring6 = ring(6.0, 0.202, 0.167, 1.0, vec2(0.002,-0.000));
ring ring7 = ring(7.0, 0.1225, 0.094, 1.0, vec2(0.002,-0.000));

ring white = ring(1.0, 0.203, 0.094+syn_BassLevel*0.015+syn_BassPresence*0.015, 1.0, vec2(0.002,-0.000));
ring pupil = ring(1.0, 0.093, 0.0, 1.0, vec2(0.002,-0.000));
ring pupilShine = ring(0.0, 0.083, 0.053, 1.0, vec2(0.002,-0.000));



float map(vec2 p) {
  for(int i = 0; i < 100; i++) {
    float a = sin(float(i * 10) + (p.x * p.y) + time) + cos(p.y);
    if(a > 0.0) return a;
  }
  return 0.0;
}


float effect( void ) {
  vec2 uv = ( gl_FragCoord.xy / resolution.xy ) * 2.0 - 1.0;
  float color = map(uv * 10.0);
  color = clamp(color-0.5, 0.0, 1.0);
  color *= 0.2;
  return color;
}

float genRing(ring ringIn, float modVal, bool modIfTrue){
  vec2 polPos = _toPolar(2.0*_uvc*vec2(ringIn.squishVal,1.0)+ringIn.trans);
  float retVal = 0.0;
  float modder = 0.0;
  if (modIfTrue){
    modder = 1.0;
  }
  vec2 polPos2 = _toPolarTrue(1.5*_uvc);
  vec2 polPos3 = _toPolarTrue(_rotate(1.5*_uvc, PI));
  vec2 mixerPos = _toPolarTrue(_rotate(1.5*_uvc, PI/2));

  float pushOut = pow((1.0+ringIn.index/7.0),2.0)*pow(syn_MidPresence,2.0)*0.05*modder;
  float ripples = pow(ringIn.radIn,1.8)*(syn_HighPresence+syn_MidHighPresence)*pow(syn_HighHits,1.7)*_rand(mix(polPos2*vec2(1.0,10.0), polPos3*vec2(1.0,3.0), abs(polPos2.y-0.5)*2.0)+syn_HighTime)*0.1;
  // pushOut = 0.0;
  if (polPos.x>ringIn.radIn+pushOut+ripples && polPos.x<ringIn.radOut+pushOut+ripples){
    retVal = 1.0;
  }
  return retVal;
}

vec3 genRingCols(vec3 ringColsIn[7], float ringSpacesIn[7], float whiteSpace){
  vec3 retCol = vec3(0.0);
  for (int i=0; i<7; i++){
    retCol += ringColsIn[i]*ringSpacesIn[i];
    whiteSpace = whiteSpace - ringSpacesIn[i];
  }
  retCol += clamp(vec3(1.0)*whiteSpace, vec3(0.0), vec3(1.0));
  return retCol;
}

float ring1Gaps[4] = float[4](0.16, 0.34, 0.585, 0.92);
float ring2Gaps[4] = float[4](0.14, 0.365, 0.57, 0.937);
float ring3Gaps[6] = float[6](0.13, 0.3, 0.57, 0.63, 0.81, 0.844);
float ring4Gaps[2] = float[2](0.645, 0.9);
float ring5Gaps[2] = float[2](0.667, 0.965);
float ring6Gaps[2] = float[2](0.175, 0.49);
float ring7Gaps[2] = float[2](0.53, 0.88);
float pupilShineGaps[2] = float[2](0.06, 0.2);

vec2 correctImageCoords(vec2 uv, vec2 imageSize){
  float wr = RENDERSIZE.x/imageSize.x;
  float hr = RENDERSIZE.y/imageSize.y;
  float ratio = _ratioFix(wr, hr);

  uv.y *= hr/ratio;
  uv.x *= wr/ratio;

  uv.x += -wr/ratio/2 + 0.5;
  uv.y += -hr/ratio/2 + 0.5;

  return uv;
}

vec4 renderLogo() {
  vec2 pos = _uvc;
  pos += vec2(0.5,0.5);
  pos = clamp(pos, vec2(0.0, 0.0), vec2(1.0, 1.0));

  vec2 polarPos = quadrantToPolar(_uvc);

  float rotTime = script_time;

  float ring1Space = genRing(ring1, 2.11, true);
  polarPos.y = fract(polarPos.y + mix(rotTime*0.02, floor(rotTime*0.02), pow(1.0-syn_FadeInOut, 3.0)));
  if (((polarPos.y > ring1Gaps[0])&&(polarPos.y < ring1Gaps[1]))||((polarPos.y > ring1Gaps[2])&&(polarPos.y < ring1Gaps[3]))) {
    ring1Space *= 0.0;
  }
  polarPos.y = fract(polarPos.y + mix(rotTime*0.02, floor(rotTime*0.02), pow(1.0-syn_FadeInOut, 3.0)));
  float ring2Space = genRing(ring2, 3.17, true);
  if (((polarPos.y > ring2Gaps[0])&&(polarPos.y < ring2Gaps[1]))||((polarPos.y > ring2Gaps[2])&&(polarPos.y < ring2Gaps[3]))) {
    ring2Space *= 0.0;
  }
  polarPos.y = fract(polarPos.y + mix(rotTime*0.02, floor(rotTime*0.02), pow(1.0-syn_FadeInOut, 3.0)));
  float ring3Space = genRing(ring3, 3.47, true);
  if (((polarPos.y > ring3Gaps[0])&&(polarPos.y < ring3Gaps[1]))||((polarPos.y > ring3Gaps[2])&&(polarPos.y < ring3Gaps[3]))||((polarPos.y > ring3Gaps[4])&&(polarPos.y < ring3Gaps[5]))) {
    ring3Space *= 0.0;
  }
  polarPos.y = fract(polarPos.y + mix(rotTime*0.04, floor(rotTime*0.04), pow(1.0-syn_FadeInOut, 3.0)));
  float ring4Space = genRing(ring4, 2.2, true);
  if ((polarPos.y > ring4Gaps[0])&&(polarPos.y < ring4Gaps[1])) {
    ring4Space *= 0.0;
  }
  polarPos.y = fract(polarPos.y + mix(rotTime*0.04, floor(rotTime*0.04), pow(1.0-syn_FadeInOut, 4.0)));
  float ring5Space = genRing(ring5, 1.3, true);
  if ((polarPos.y > ring5Gaps[0])&&(polarPos.y < ring5Gaps[1])) {
    ring5Space *= 0.0;
  }
  polarPos.y = fract(polarPos.y + mix(rotTime*0.08, floor(rotTime*0.08), pow(1.0-syn_FadeInOut, 4.0)));
  float ring6Space = genRing(ring6, 1.1, true);
  if ((polarPos.y > ring6Gaps[0])&&(polarPos.y < ring6Gaps[1])) {
    ring6Space *= 0.0;
  }
  polarPos.y = fract(polarPos.y + mix(rotTime*0.16, floor(rotTime*0.16), pow(1.0-syn_FadeInOut, 4.0)));
  float ring7Space = genRing(ring7, 0.537, true);
  if ((polarPos.y > ring7Gaps[0])&&(polarPos.y < ring7Gaps[1])) {
    ring7Space *= 0.0;
  }

  float whiteSpace = genRing(white, 0.0, false);
  float pupilSpace = genRing(pupil, 0.0, false);
  float pupilShineSpace = genRing(pupilShine, syn_HighHits+1.0, false);
  if ((polarPos.y < pupilShineGaps[0])||(polarPos.y > pupilShineGaps[1])) {
    pupilShineSpace *= 0.0;
  }
  float ringSpaces[7] = float[7](ring1Space, ring2Space, ring3Space, ring4Space, ring5Space, ring6Space, ring7Space);

  vec3 ringCol = genRingCols(colRings, ringSpaces, whiteSpace);
  ringCol += pupilShineSpace;

  vec3 colLogo = texture(logo, clamp(pos*2.0-vec2(0.5),0.0,1.0)).rgb;
  // colLogo *= (0.5);


  vec3 finalCol = mix(ringCol, colLogo, no_motion);

  return vec4(finalCol, 1.0);
}

vec4 renderMain(){
  colRings = vec3[7](_normalizeRGB(46, 191, 239), _normalizeRGB(97, 164, 147), _normalizeRGB(102, 84, 164), 
  _normalizeRGB(238, 38, 141), _normalizeRGB(245, 135, 56), _normalizeRGB(248, 239, 36), _normalizeRGB(46, 191, 239));

  if (PASSINDEX == 0.0){
    return renderLogo();
  } else if (PASSINDEX == 1.0){
    vec2 coords = correctImageCoords(_uv*4.5-vec2(1.44,0.55)*(4.5/3.7), textureSize(text, 0));
    vec3 textCol = texture(text, clamp(coords,0.0,1.0)).rgb;
    if ((coords.x > 0.91)||(coords.x < 0.09)||(coords.t < 0.2)||(coords.t > 0.8)){
      textCol = vec3(0.0);
    }
    return texture(logoPass, _uv-vec2(0.0,0.15))+vec4(textCol,0.0)*mix(1.0, 0.0, pow(syn_FadeInOut,0.2)*(1.0-no_motion));
  }
}
