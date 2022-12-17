#define MAX_ITER 75
float noise( in vec2 x )
{
	return sin(1.25*x.x)*sin(1.25*x.y);
}
mat2 m = mat2(.85,.65,-.65,.85);
vec2 circNoise = vec2((TIME), (TIME));

float fbm6( vec2 p )
{
    p += pow(syn_BassLevel*0.75, 2.);
    float f = 0.0;
    f += 0.500000*(0.5+0.5*noise( p )); p = m*p*2.02;
    f += 0.250000*(0.5+0.5*noise( p )); p = m*p*2.03;
    f += 0.125000*(0.5+0.5*noise( p )); p = m*p*2.01;
    f += 0.062500*(0.5+0.5*noise( p )); p = m*p*2.04;
    f += 0.031250*(0.5+0.5*noise( p )); p = m*p*2.01;
    f += 0.015625*(0.5+0.5*noise( p ));
    return f/(0.96875);
}

float cameraZ = fly_in_out + smoothTime*1.50;

float lengthRep = 1500;

vec3 palette( in float t)
{
vec3 a = vec3(0.600, 0.500, 0.500);
vec3 b = vec3(1.000, 0.500, 0.500);
vec3 c = vec3(0.500, 0.690, 0.750);
vec3 d = vec3(0.500, 0.500, 0.500);
  
    return a + b*cos( 6.28318*(c*t+d) );
}

vec3 filter_() {
  vec2 delta = 1. / RENDERSIZE;
  
  vec3 val = texture(backbuffer, _uv).xyz;


  vec3 l = texture(backbuffer, _uv + vec2(0., delta.y)).xyz;
  vec3 r = texture(backbuffer, _uv - vec2(0., delta.y)).xyz;
  vec3 u = texture(backbuffer, _uv + vec2(delta.x, 0.)).xyz;
  vec3 d = texture(backbuffer, _uv - vec2(delta.x, 0.)).xyz;

  vec3 n = vec3(_rand(_uvc+fract(TIME))) - 0.5;
  
  vec3 bloom = max(val, max(max(l, r), max( u, d)));
  bloom = bloom  + l + r + u + d;
  bloom /= 5.; // orlando;
  return bloom + n/9.;

}

vec3 rotateX(vec3 p, float angle)
{
    float c = cos(angle);
    float s = sin(angle);
    return vec3(p.x, c*p.y+s*p.z, -s*p.y+c*p.z);
}

vec3 rotateY(vec3 p, float angle)
{
    float c = cos(angle);
    float s = sin(angle);
    return vec3(c*p.x-s*p.z, p.y, s*p.x+c*p.z);
}

vec3 rotateZ(vec3 p, float angle)
{
    float c = cos(angle);
    float s = sin(angle);
    return vec3(c*p.x+s*p.y, -s*p.x+c*p.y, p.z);
}

mat3 rotationMatrix(vec3 axis, float angle) {
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}

vec3 opRep( vec3 p, vec3 c )
{
    vec3 q = mod(p, c) - 0.5 * c;
    
    return q;
}

float sdCappedCylinder( vec3 p, vec2 h )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float opS( float d1, float d2 )
{
    return max(-d1,d2);
}

float sdCylinder( vec3 p, vec3 c )
{
  return length(p.xz-c.xy)-c.z;
}

float hex(vec3 p, vec2 h) {
  vec3 q = abs(p);
  return max(q.z-h.y, max((q.x*0.866025+q.y*0.5),q.y) - h.x);
}

float sdTriPrism( vec3 p, vec2 h )
{
    vec3 q = abs(p);
    return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
}
vec2 polar(vec2 c) {
  return vec2(atan(c.y, c.x), length(c));
}
vec2 cart(vec2 p) {
  return p.y * vec2(cos(p.x), sin(p.x));
}

// polynomial smooth min (k = 0.1);
float smin( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float sdSphere( vec3 p, float s )
{
  return length(p)-s;
}

float impulse( float k, float x )
{
    float h = k*x;
    return h*exp(1.0-h);
}  

float cubicPulse( float c, float w, float x )
{
    x = abs(x - c);
    if( x>w ) return 0.0;
    x /= w;
    return 1.0 - x*x*(3.0-2.0*x);
}

vec3 DE(vec3 pos) {
  float dist = 10000.0;
  float obj = 0.0;

  vec3 staticPos = vec3(pos.xy, pos.z-cameraZ);
  staticPos.xy +=_uvc;
  pos.xy+=_uvc;
  pos.y -= 10.0*_pulse(staticPos.z*0.1, (1.0-wave_pulse)*10.0, 4.0)*sqrt(wave_pulse);

  //Corridor
  vec3 corridorPos = pos;
  corridorPos.xy = _rotate(corridorPos.xy+_uvc*FOV, corridorPos.z*2*PI/30.0);
  float corridor = 15.0 - abs(corridorPos.y);

  if (undulating > 0.5){
    corridor -= (0.5+0.4*sin(pos.z*0.05+0.05*smoothTimeC*PI))*10.0;
  }

  dist = min(corridor, dist);

  // //COLUMNS!
  vec3 columnPos = pos;

  float columns = 1000.0;
  vec3 columnPos2 = pos;

  if (many_beams>0.5){
    columnPos = rotateZ(pos, 0.1);
    columnPos.y += cubicPulse(30.0+20.0*sin(syn_BassTime*0.25), 10.0, staticPos.z)*4.0*syn_BassLevel;
    columnPos = opRep(columnPos, vec3(2.5-syn_HighPresence*1.0,10.0,20.0));
    columnPos = rotateZ(columnPos, columnPos.z*0.1);
    columns = sdTriPrism(columnPos, vec2(0.1+0.4*syn_HighHits, 100.0));
  } else {
    columnPos2 = rotateZ(columnPos2, 0.1);
    columnPos2 = opRep(columnPos2, vec3(5.0,20.0,20.0));
    columns = sdTriPrism(columnPos2+vec3(0.0,0.0,sin(pos.z*0.1)*10.0), vec2(0.0, 100.0));
    // columns = mix(1000, columns, columns_on);

  }
    // dist = smin(dist, columns, 5.0);
  // } else {
  //   columns = mix(1000, columns, columns_on);
    // columns = mix(1.0, columns, columns_on);
    dist = mix(smin(dist, columns, 1.5+syn_Presence*3.5), dist, 1.0-beams_on);



  // }

  ////Sphere
  // vec3 sphereLightPos = staticPos + vec3(sin(TIME),cos(TIME),-5.0);

  // float sphereLight = sdSphere(sphereLightPos, 0.5+syn_HighHits);
  // dist = smin(dist, sphereLight, 3.0);


  //CUTS!
  //spherePos = opRep(spherePos, vec3(1.0,1.0,1.0));
  //spherePos += vec3(sin(TIME), cos(TIME), 0.0);
  // float sphere = sdTriPrism(spherePos, vec2(1.0));
  // dist = min(corridor, sphere);

  obj = clamp(1.0-columns,0.0,1.0);
  
  float glow = 1/columns;

  return vec3(dist, obj, glow);
}

vec4 bumpMapParams1 = vec4(2.0,7.0,0.01,-0.01);
vec4 bumpMapParams2 = vec4(2.0,3.0,-0.01,0.01);

// float bumpmap(vec2 uv)
// {
//     float b1 = bumpMapParams1.x*(1.0 - texture(iChannel0, (uv + iGlobalTime*bumpMapParams1.zw)*bumpMapParams1.y).x);
//     float b2 = bumpMapParams2.x*(1.0-texture(iChannel0, (uv + iGlobalTime*bumpMapParams2.zw)*bumpMapParams2.x).x);
//     return b1+b2;
// }


float bumpmap(vec2 uv)
{
  return sin(uv.x*100.0+uv.y*500.0);
    // float b1 = bumpMapParams1.x*(1.0 - texture(bump_texture, (uv)*bumpMapParams1.y - vec2(syn_HighTime*0.5,0.0)).x);
    // float b2 = bumpMapParams2.x*(1.0 - texture(bump_texture, (uv)*bumpMapParams2.x - vec2(syn_HighTime*0.5,0.0)).x);
    // return b1+b2;
}

vec3 doBumpMap(vec2 uv, vec3 nor, float bumpfactor)
{
   
    const float eps = 0.001;
    float ref = bumpmap(uv); 
    
    vec3 grad = vec3(bumpmap(vec2(uv.x-eps, uv.y))-ref, 0.0, bumpmap(vec2(uv.x, uv.y-eps))-ref); 
             
    grad -= nor*dot(nor, grad);          
                      
    return normalize( nor + grad*bumpfactor );
}

vec3 gradient(vec3 p, float t) {
      vec2 e = vec2(0., t);

      return normalize( 
        vec3(
          DE(p+e.yxx).x - DE(p-e.yxx).x,
          DE(p+e.xyx).x - DE(p-e.xyx).x,
          DE(p+e.xxy).x - DE(p-e.xxy).x
        )
      );
    }

vec3 raycast() {
  
  vec3 camera = vec3( 0., 0., cameraZ);
  camera.z = mod(camera.z, lengthRep);

  vec2 screenPos = -1.0 + 2.0 * _uv;
  screenPos.x *= RENDERSIZE.x / RENDERSIZE.y;
  screenPos.xy += FOV*((_uvc*PI)+camera.xy);

  screenPos.xy *= 1.0+WarpMirror*((_uvc*PI*PI+_uv)-1.);

  vec3 ray = normalize(vec3( screenPos, 1.0));
  // ray = rotateX(ray, TIME);

  ray = rotateX(ray, sin(smoothTime*PI*0.01225)*0.125);
  ray = rotateY(ray, sin(smoothTime*PI*0.0125)*0.125);
  //ray = rotateZ(ray, smoothTime*0.025*rotating);
    ray += 0.5*(fbm6( (ray.xy*_uvc))*pow(basshits*0.5, 1.4))*Distortion;

  ray.xz =  _rotate(ray.xz, LookXY.x*PI*PI+_uvc.y*PI*Perspective.x);
  ray.yz = _rotate(ray.yz, LookXY.y*PI+_uvc.y*PI*Perspective.y);
  ray.xy = _rotate(ray.xy, Rotate*PI);

  // raycasting parameter 
  float t  = 0.;
  vec3 pos;
  int iter = 0;
  bool hit = false;
  vec3 dist;
  float glow = 0.0;
  float minDist = pow(0.6+0.0*5.0,10.0);

  // ray stepping 
  for(int i = 0; i < MAX_ITER; i++) {
    pos = camera + ray * t;
    dist = DE(pos);
    glow += dist.z;

    if (abs(dist.x) < minDist) {
      hit = true;
      break;
    }
              
    t += dist.x * 0.5;//step_size
    iter ++;

  }

  float obj = dist.y;
  vec3 norm = gradient(pos, 0.001);
  vec2 uv = _rotate(pos.xy, pos.z*2*PI/30.0);

  vec3 bump = doBumpMap(vec2(uv.x*0.01, pos.z*0.01), norm, 4.0);
  bump = mix(bump, vec3(0.0), clamp(obj,0.0,1.0))*(1.0+0.25*sin(smoothTimeC*0.0125));
  bump = mix(vec3(0.0), bump, walls_texture);



  vec3 topLightPos = -ray + vec3(0.0,-1.0,-2.0);
  float topLight = clamp(dot(normalize(norm+bump), topLightPos), 0.0, 1.0);

  vec3 botLightPos = -ray + vec3(sin(smoothTimeB*0.25),cos(smoothTimeB*0.25),2.0+sin(smoothTimeB*0.25));
  float botLight = clamp(dot(normalize(norm+bump), botLightPos), 0.0, 1.0);

  float generalLight = clamp(dot(normalize(norm+bump), -ray), 0.0, 1.0);

    // float generalLight = clamp(dot(normalize(norm+bump), -ray), 0.0, 1.0);


  vec3 lp = camera + vec3(0.0,0.0,0.5+0.5*sin(script_time*PI*0.25)*200.0); // Light position - Back from the screen.
  // vec3 lp = camera; // Light position - Back from the screen.

  vec3 ld = lp - pos;
  float lDist = max(length(ld), 0.001);
  ld /= lDist;

  float spec = pow(max(dot(reflect(-ld, normalize(norm+bump)), -ray), 0.), 24.); 

  float ambientOcclusion = 1. -  float(iter) / float(MAX_ITER);
  
  vec3 color = vec3(0.);

  //vec3 lightCol = vec3(0.1,0.0,0.0);
  vec3 lightCol = vec3(0.,0.0,0.0);
  
  // lightCol = cross(lightCol, pos);
  // lightCol = vec3(1.0);

  if (lights_react > 0.5){
    topLight *= (highhits);
    botLight *= (syn_MidLevel*0.5+syn_MidHighLevel*0.5);
  }

  vec3 columnsCol = vec3(0.0);
  vec3 cavernCol = vec3(0.0);

  vec3 surfaceCol = vec3(1.0);

  if (colorful>0.5){
  //  surfaceCol = _palette(dot(normalize(norm+bump), -ray), vec3(1.020, 0.520, 0.500), vec3(0.500, 0.500, 0.500), vec3(0.500, 0.240, 0.750), vec3(0.480, 0.090, 0.870));
    surfaceCol = _palette(dot(normalize(norm+bump), -ray), vec3(1.020, 0.520, 0.500), vec3(0.500, 0.500, 0.500), vec3(0.500, 0.240, 0.750), vec3(0.480, 0.1090, 0.870));
  
  }

  if (hit){
    columnsCol += lightCol*botLight;
    cavernCol += spec*vec3(1.0);

    color = mix(vec3(0.0), columnsCol, clamp(obj,0.0,1.0));
    color = mix(color, cavernCol, clamp(obj, 1.0, 2.0)-1.0);
    if (light_chaser>0.5){
      color += generalLight*vec3(1.0)*_pulse(sin((pos.z-cameraZ)*0.01), fract(smoothTime/1.0), 0.1);
    }
    color += botLight*surfaceCol;

    color += spec*surfaceCol;
    color *= pow(ambientOcclusion,0.1);
  }

  float fogAmount = 1.0 - exp( -t*0.005);
  vec3  fogColor  = vec3(0.3+syn_BassPresence*0.6,0.6,0.7);
  color = mix(color, fogColor, fogAmount*syn_Presence);

  color += clamp(pow(glow*beam_glow,0.5)*0.1,0.0,1.0);

  color = pow(color, vec3(1.5));

  return color;


}

vec4 renderMain() {
  vec4 out_FragColor = vec4(0.0);

  if (PASSINDEX == 0) {
    out_FragColor = vec4(raycast(), 1.0);
  } else if (PASSINDEX == 1) {
    // out_FragColor = texture(backbuffer, _uv);
    out_FragColor = vec4(filter_(), 1.0);

    vec4 bCol = vec4(0.4, 0.7, 0.9, 1.0);
    vec4 oCol = vec4(0.9, 0.8, 0.6, 1.0);

    out_FragColor = mix(out_FragColor, pow(out_FragColor,vec4(2.0))*bCol, clamp(abs(_uvc.y)-0.3, 0.0, 1.0));
    out_FragColor = mix(out_FragColor, pow(out_FragColor,vec4(2.0))*oCol, clamp(0.7-abs(_uvc.y-0.3)*2.0, 0.0, 1.0));

  }

  return out_FragColor;
}
