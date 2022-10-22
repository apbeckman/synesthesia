

// Alien engine
// Based upon: https://www.shadertoy.com/view/4ds3zn

#define PI         3.141592654
#define TAU        (2.0*PI) 
#define TOLERANCE  0.0003
#define REPS       11
#define MAX_DIST   25.
#define MAX_ITER   120

const vec3  green  = vec3(1.75, 2.0, 1.40);
const vec3  dark   = vec3(0.2);

void rot(inout vec2 p, float a) {
  float c = cos(a);
  float s = sin(a);
  p = vec2(c*p.x + s*p.y, -s*p.x + c*p.y);
}

void r45(inout vec2 p) {
	p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

float apollian(vec3 p, float tolerance, out int layer) {
  const float s = 1.9;
  float scale = 1.0;

  float r = 0.2;
  vec3 o = vec3(0.22, 0.0, 0.0);

  float d = 0.0;

  for(int i = 0; i < REPS; ++i) {
    p = (-1.00+ 2.0*fract((0.5+test)*p+0.5));
//    rot(p.xz, -float(i)*PI/4.0);
    r45(p.xz);

    float r2 = dot(p,p) + 0.0;
    float k = s/r2;
    float ss = pow((1.0 + float(i)), -0.15+sin(smoothTime*0.1)*0.05-.075);
    p *= pow(k, ss);
    scale *= pow(k, -ss*ss);
    d = (0.25)*abs(p.y)*scale;
    layer = i;
    if(abs(d) < tolerance) break;
  }
  
  return d;
}

float df(vec3 p, float tolerance, out int layer) { 
  float d = apollian(p, tolerance, layer);
  return d;
} 


float intersect(vec3 ro, vec3 rd, out int iter, out int layer) {
  float res;
  float t = 1.6;
  iter = MAX_ITER;
    
  for(int i = 0; i < MAX_ITER; ++i) {
    vec3 p = ro + rd * t;
    float tolerance = TOLERANCE * t;
    res = df(p, tolerance, layer);
    if(res < tolerance || res > MAX_DIST) {
      iter = i;
      break;
    }
    t += res;
  }
    
  if(res > MAX_DIST) t = -1.;
    
  return t;
}

float ambientOcclusion(vec3 p, vec3 n) {
  float stepSize = 0.01;
  float t = stepSize;

  float oc = 0.0;
  
  int layer;

  for(int i = 0; i < 12; i++) {
    float tolerance = TOLERANCE * t;
    float d = df(p + n * t, tolerance, layer);
    oc += t - d;
    t += stepSize;
  }

  return clamp(oc, 0.0, 1.0);
}

vec3 normal(in vec3 pos) {
  vec3 eps = vec3(.001,0.0,0.0);
  vec3 nor;
  int layer;
  float tolerance = TOLERANCE * eps.x;
  nor.x = df(pos+eps.xyy, tolerance, layer) - df(pos-eps.xyy, tolerance, layer);
  nor.y = df(pos+eps.yxy, tolerance, layer) - df(pos-eps.yxy, tolerance, layer);
  nor.z = df(pos+eps.yyx, tolerance, layer) - df(pos-eps.yyx, tolerance, layer);
  return normalize(nor);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
 
  vec2 q=fragCoord.xy/RENDERSIZE.xy; 
  vec2 uv = -1.0 + 2.0*q; 
  uv.x*=RENDERSIZE.x/RENDERSIZE.y; 
    
  vec3 la = vec3(0.0,0.5,0.0); 
  vec3 ro = vec3(2.5, 1.5, 0.0);
  rot(ro.xz, smoothTime/20.0);

  vec3 cf = normalize(la-ro); 
  vec3 cs = normalize(cross(cf,vec3(0.0,1.0,0.0))); 
  vec3 cu = normalize(cross(cs,cf)); 
  vec3 rd = normalize(uv.x*cs + uv.y*cu + 3.0*cf);

  vec3 bg = mix(dark*0.25, dark*0.5, smoothstep(-1.0, 1.0, uv.y));
  vec3 col = bg;

  vec3 p=ro; 

  int iter = 0;
  int layer = 0;
  
  float t = intersect(ro, rd, iter, layer);
    
  if(t > -1.0) {
    p = ro + t * rd;
    vec3 n = normal(p);
    float fake = float(iter)/float(MAX_ITER);
    float fakeAmb = exp(-fake*fake*4.0);
    float amb = ambientOcclusion(p, n);
    

    vec3 dif;

    float ll = length(p);

    if (layer == 0)
    {
      dif = 0.75*green;
      //dif += highhits;

    } else {
      dif = green*pow((1.0 + 0.5*cos(-PI*(2.0+highhits*0.5)*float(layer)/float(REPS) + smoothTimeB*0.5 - 0.5*PI*(ll))), 4.0)/pow(float(layer), 1.5);

    }


    const float fogPeriod = TAU*2.0;
    float fogHeight = 0.25 + 0.325*(abs(p.y) + 0.125*(sin(fogPeriod*p.x) * cos(fogPeriod*p.z)));
    float dfog = (fogHeight - ro.y)/rd.y;
    float fogDepth = t > dfog && dfog > 0.0 ? t - dfog : 0.0;
    float fogFactor = exp(-fogDepth*4.0);
  fogFactor -= highhits*0.125;
    col = dif;
    col *= vec3(mix(1.0, 0.125, pow(amb, 3.0)))*vec3(fakeAmb);
    col = mix(green*0.5, col, fogFactor); 
    col = mix(bg, col, exp(-0.0125*t*t)); 
  } 

  //float pp = 1.0 - (1.0 - step(0.5, q.y))*smoothstep(0.85, 1.3, length(2.0*q-1.0));
  
    
  col;
    

  fragColor=vec4(col.x,col.y,col.z,1.0); 
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}