

#define iResolution vec3(RENDERSIZE, 1.)
#define iTime TIME



float time;
float  timeB = smoothTimeB * 1.0 + 137.0;
float  timeC = smoothTimeC * 1.0 + 137.0;

mat2 rot(float a) {
  float ca=cos(a);
  float sa=sin(a);
  return mat2(ca,sa,-sa,ca);  
}

float box(vec3 p, vec3 s) {
  p=abs(p)-s;
  return max(p.x, max(p.y,p.z));
}

vec3 fr(vec3 p, float t) {

  //float s = 1.0 - exp(-fract(time*1.0))*0.8;
  float s = 0.7 - smoothstep(0.0,1.0,abs(fract(timeC*0.1)-0.5)*2.0)*0.3;
  for(int i=0; i<7; ++i) {
    
    float t2=t+float(i);
    p.xy *= rot(t2);
    p.yz *= rot(t2*.7);
    
    float dist = 10.0;
    p=(fract(p/dist-.5)-.5)*dist*BILLOW2;
    p=abs(p);
    p-=s;
    
  }
  
  return p;
}

float at = 0.;
float at2 = 0.;
float at3 = 0.;
float map(vec3 p) {
  
  
  vec3 bp=p;
  
  p.xy *= rot((p.z*0.023+timeC*0.1)*0.3);
  p.yz *= rot((p.x*0.087)*0.4);
  
  float t=time*0.5;
  vec3 p2 = fr(p, t * 0.2);
  vec3 p3 = fr(p+vec3(5,0,0), t * 0.23);
  
  float d1 = box(p2, vec3(1,1.3,4));
  float d2 = box(p3, vec3(3,0.7,0.4));
  
  float d = max(abs(d1), abs(d2))-0.2;
  float dist = 1.;
  vec3 p4=(fract(p2/dist-.5)-.5)*dist;
  float d3 = box(p4, vec3(0.4));
  //d = max(d, -d3);
  d = d - d3*0.4;
  
  //d = max(d, length(bp)-15);
  
  
  //float f=p.z + time*4;
  //p.x += sin(f*0.05)*6;
  //p.y += sin(f*0.12)*4;
  //d = max(d, -length(p.xy)+10);
  
  at += 0.13/(0.13+abs(d));
  
  float d5 = box(bp, vec3(4));
  
  float dist2 = 8.;
  vec3 p5=bp;
  p5.z = abs(p5.z)-13.;
  p5.x=(fract(p5.x/dist2-.5)-.5)*dist2;
  float d6 = length(p5.xz)-1.;
  
  at2 += 0.2/(0.15+abs(d5));
  at3 += 0.2/(0.5+abs(d6));
  
  return d;
}

void cam(inout vec3 p) {
  
  float t=timeC*0.05;
  p.yz *= rot(t);
  p.zx *= rot(t*1.2);
}

float rnd(vec2 uv) {  
  return fract(dot(sin(uv*752.322+uv.yx*653.842),vec2(254.652)));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{    
    
  time = smoothTime * 1.0 + 137.0;
    
  vec2 uv = vec2(fragCoord.x / iResolution.x, fragCoord.y / iResolution.y);
  uv -= 0.5;
  uv /= vec2(iResolution.y / iResolution.x, 1);
  
  float factor = 0.9 + 0.1*rnd(uv);
  //factor = 1;

  vec3 s=vec3(0,0,-15);
  vec3 r=normalize(vec3(-uv, 1));
  
  cam(s);
  cam(r);
  
  vec3 p=s;
  int i=0;
  
  for(int i=0; i<70; ++i) {
    float d=abs(map(p));
    d = abs(max(d, -length(p-s)+6.));
    d *= factor * BILLOW;
    if(d<0.001) {
      d = 0.1;
      //break;
    }
    p+=r*d;
  }
  
  vec3 col=vec3(0);
  //col += pow(1-i/101.0,8);
  
  vec3 sky = mix(vec3(1,0.5,0.3), vec3(0.2,1.5,0.7), pow(abs(r.z),8.));
  sky = mix(sky, vec3(0.4,0.5,1.7), pow(abs(r.y),8.));
  
  //col += at*0.002 * sky;
  col += pow(at2*0.008, 1.) * sky;
  col += pow(at3*0.072, 2.) * sky * vec3(0.7,0.3,1.0) * 2.;
  
  col *= 1.2-length(uv);
  
  col = 1.0-exp(-col*15.0);
  col = pow(col, vec3(1.2));
  col *= 1.2;
  //col += 0.2*sky;
  
  //col = vec3(rnd(uv));
  
  fragColor = vec4(col, 1);
}

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);

    mainImage(out_FragColor, _xy.xy);

return out_FragColor; 
 } 
