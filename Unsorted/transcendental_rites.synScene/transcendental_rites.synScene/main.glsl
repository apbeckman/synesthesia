//Transcendental Rite - Psybernautics - Alex Tiemann - 2021
#define FAR 80.
// Tri-Planar blending function. Based on an old Nvidia writeup:
// GPU Gems 3 - Ryan Geiss: https://developer.nvidia.com/gpugems/GPUGems3/gpugems3_ch01.html

vec3 tex3D(sampler2D t, in vec3 p, in vec3 n ){

    

    n = max(abs(n) - .2, 0.001);

    n /= dot(n, vec3(1));

    vec3 mirP = p;

    vec3 modP = mod(p, vec3(2.0));

    if (modP.x > 1.0){

        mirP.x = 1.0-p.x;

    }

    if (modP.y > 1.0){

        mirP.y = 1.0-p.y;

    }

    if (modP.z > 1.0){

        mirP.z = 1.0-p.z;

    }

    vec3 tx = _contrast(_invertImage(vec4(texture(t, mirP.yz).xyz, 0.0)), _Media_Contrast).rgb;

    vec3 ty = _contrast(_invertImage(vec4(texture(t, mirP.zx).xyz, 0.0)), _Media_Contrast).rgb;

    vec3 tz = _contrast(_invertImage(vec4(texture(t, mirP.xy).xyz, 0.0)), _Media_Contrast).rgb;


    return (tx*tx*n.x + ty*ty*n.y + tz*tz*n.z);

}

float smin( float a, float b, float k )
{
    float h = max(k-abs(a-b),0.0);
    return min(a, b) - h*h*0.25/k;
}

vec2 smin( vec2 a, vec2 b, float k )
{
    float h = clamp( 0.5+0.5*(b.x-a.x)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

mat2 r2d(float a) 
{
    float c = cos(a), s = sin(a);
    return mat2(c, s, -s, c);
}

vec3 re(vec3 p, float d) 
{return mod(p - d * .5, d) - d * .5;}

void amod(inout vec2 p, float d) 
{
    float a = re(vec3(atan(p.y, p.x)), d).x;
    p = vec2(cos(a), sin(a)) * length(p);
}

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

vec2 map(vec3 p) {

    vec2 scene = vec2(0.0);
    vec3 u = p;
    u.y -= 1.0;
    u = p;
    u.xz *= r2d(geo_dynamic_time*0.0-u.y*0.0625);
    u.xz += normalize(u.xz)*(0.04)*(sin( u.y*0.01325));
    u.x -= 2.0*sin(u.y*0.125+pBass);
    u.x += sin(u.y*0.0625 + pMid);
    amod(u.xz, 2.0*PI /  floor(mix(spire_count, 2.0 + syn_RandomOnBeat*4.0, rand_spire_count )));
    u.x -= 15.0+1.5*sin(u.y*0.0625+geo_dynamic_time);
    u.x += sin(u.y*0.0625 + geo_dynamic_time);
    u.xz *= r2d(u.y*0.5);
    amod(u.xz, 2.0*PI / clamp(floor(spire_count) / 2.0, 3.0, 16.0));
    u.x -= 3.0+ sin(u.y*0.125+geo_dynamic_time*0.75);
    u.x += sin(u.y*0.0625 + geo_dynamic_time);
    float toobs = length(u.xz) - 1.0;
    scene.x = toobs;
    return scene;
}


vec2 mapRef(vec3 p) {


    vec2 scene = vec2(2.0);

    scene.x = length(p*0.125) - 2.0;
    return scene;
    
}

vec2 march(vec3 ro, vec3 rd){

    vec2 t = vec2(0.0);
    vec2 d = vec2(0.0);

    for (int i = 0; i < 120; i++){

        d = map(ro + rd*t.x);
        if(abs(d.x)<.0005 || t.x>FAR) break;        
        t.x += d.x*.75;
    }

    t.y = d.y;
    return t;

}

vec2 marchRef(vec3 ro, vec3 rd){

    vec2 t = vec2(0.0);
    vec2 d = vec2(0.0);

    for (int i = 0; i < 8; i++){

        d = mapRef(ro + rd*t.x);

        if(abs(d.x)<.002 || t.x>(FAR*1.5)) break;
        t += d;

    }
    return t;

}


// Tetrahedral normal by IQ.

vec3 getNormal( in vec3 p ){

    vec2 e = vec2(.0025, -.0025); 

    return normalize(

        e.xyy * map(p + e.xyy).x + 

        e.yyx * map(p + e.yyx).x + 

        e.yxy * map(p + e.yxy).x + 

        e.xxx * map(p + e.xxx).x);

}



vec3 getObjectColor(vec3 p, vec3 n, float mat){

    vec3 col = base_color;
    col = mix(col,_hsv2rgb(vec3((sdBox(p*0.01+n.y*0.25, vec3(dynamic_time*0.25))), 1.0, 1.0)),rainbow);
    col = syn_MediaType > 0.5 ? mix(col, mix(_loadUserImage().xyz, tex3D(syn_UserImage, p*0.03125, n), media_uv_or_3d)*1.5, media_color) : col;

    return col;

}


vec3 doColor(in vec3 sp, in vec3 rd, in vec3 sn, in vec3 lp, vec2 t){

    vec3 ld = lp-sp;
    float lDist = max(length(ld), .001);
    ld /= lDist;
    float atten = 2.0 / (1. + lDist*.2 + lDist*lDist*.1);
    float diff = max(dot(sn, ld), 0.);
    float spec = pow(max( dot( reflect(-ld, sn), -rd ), 0.), 2.0);
    vec3 objCol = getObjectColor(sp, sn, t.y);
    float lightStrength = specular_strength;
    vec3 sceneCol = (objCol*(diff*lightStrength*0.25+30.0 + (.15+media_color*5.5)) + (vec3(1.0)-objCol)*spec*lightStrength) * atten;
    float fogF = smoothstep(0., .95, t.x/FAR);
    sceneCol = mix(sceneCol, vec3(0.0), fogF); 

    return sceneCol;
}

mat3 cam(vec3 ro, vec3 ta)
{
    vec3 f=normalize(ta-ro);
    vec3 r=normalize(cross(f,vec3(0.,1.,0.)));
    vec3 u=normalize(cross(r,f));
    return mat3(r,u,f);
}

vec4 renderMainImage() {

    vec3 ro = vec3(0.0,0.0,cos(dynamic_time));
    vec3 lk = vec3(0.0,5.0,0.0);
    ro += vec3(0.0, 0.0, -40.0);
    vec3 lp = ro;
    lp.y -= inTube*0.35;
    ro.z += inTube;
    lp.z = 10.0;
    
    ro.xz = _rotate(ro.xz, 0.5*dynamic_time);
    lp.xz = _rotate(lp.xz, 0.5*dynamic_time);

    float f_o_v = inTube > 38.0 ? clamp(FOV, 0.03, 0.1) : clamp(FOV, 0.5, 3.0);
	vec3 rd=cam(ro,lk)*normalize(vec3(_uvc, f_o_v ));

    //Based on Perspective Camera - circuit bending - meebs
    float T = script_time;
	vec2 v = (_xy.xy / RENDERSIZE.xy) + RENDERSIZE.y;
	v.y -= 0.5;
    float pi = PI;
    float twpi = 2.0*PI;
	float th =  v.y * pi, ph = v.x * twpi;
    vec3 ssp = vec3( sin(ph) * cos(th), sin(th), cos(ph) * cos(th) );
    float time = dynamic_time;
    ssp.xy = _rotate(ssp.xy, pMid);
    ssp.yz = _rotate(ssp.yz, pMidHigh);
    ssp.xz = _rotate(ssp.xz, pBass);
    rd = mix(rd, normalize(ssp), cam_morph);

    vec2 t = march(ro, rd);
    vec3 sp = ro + rd*t.x;
    ro += rd*t.x;   
    vec3 sn = getNormal(ro);
    rd = reflect(rd, sn);
    vec2 tSave = t;
    t = marchRef(ro +  rd*.01, rd);
    ro += rd*t.x;
    sn = getNormal(ro);

    vec3 sceneColor = doColor(ro, rd, sn, lp, tSave);

    //vec3 uv = n;
    vec2 uv = _uv;
    uv.y -= _noise(uv.x*20.0 + beatStripe)*0.05 + beatStripe;
    float stripe = fract(uv.y*4.0);
    sceneColor = mix(sceneColor, sceneColor*stripe + vec3(0.0)*(1.0 - stripe), beat_stripe);

    vec2 u = _xy/RENDERSIZE.xy;
    sceneColor = mix(vec3(0), sceneColor, pow( 16.0*u.x*u.y*(1.0-u.x)*(1.0-u.y) , .5)*.5 + .5);
    vec4 color = vec4(sqrt(clamp(sceneColor, 0., 1.)), 1);
    color = pow(color, vec4(2.0));

    return color;
 } 

/**

 * Detects edges using the Sobel equation

 * @name syn_pass_edgeDetectSobel

 * @param  {sampler2D} smp texture you wish to detect edges on

 * @returns {float} edges

 */



vec4 sobelIntensity(in vec4 color){

  return color;

  //return vec4(sqrt((color.x*color.x)+(color.y*color.y)+(color.z*color.z)));

}

vec4 sobelHelper(float stepx, float stepy, vec2 center, sampler2D tex){

  // get samples around pixel

  vec4 tleft = sobelIntensity(texture(tex,clamp(center + vec2(-stepx,stepy), 0.0, 1.0)));

  vec4 left = sobelIntensity(texture(tex,clamp(center + vec2(-stepx,0), 0.0, 1.0)));

  vec4 bleft = sobelIntensity(texture(tex,clamp(center + vec2(-stepx,-stepy), 0.0, 1.0)));

  vec4 top = sobelIntensity(texture(tex,clamp(center + vec2(0,stepy), 0.0, 1.0)));

  vec4 bottom = sobelIntensity(texture(tex,clamp(center + vec2(0,-stepy), 0.0, 1.0)));

  vec4 tright = sobelIntensity(texture(tex,clamp(center + vec2(stepx,stepy), 0.0, 1.0)));

  vec4 right = sobelIntensity(texture(tex,clamp(center + vec2(stepx,0), 0.0, 1.0)));

  vec4 bright = sobelIntensity(texture(tex,clamp(center + vec2(stepx,-stepy), 0.0, 1.0)));



  vec4 x = tleft + 2.0*left + bleft - tright - 2.0*right - bright;

  vec4 y = -tleft - 2.0*top - tright + bleft + 2.0 * bottom + bright;

  vec4 color = sqrt((x*x) + (y*y));

  return color;

}

vec4 edgeDetectSobel(sampler2D tex){

	float stepSize = 1.0+syn_BassLevel*2.0;

  vec2 uv = _uv;

  if (uv.x < 0.0 || uv.y < 0.0 || uv.y > 1.0 || uv.x > 1.0) {

    return vec4(0.0);

  } 



  return sobelHelper(stepSize/RENDERSIZE.x, stepSize/RENDERSIZE.y, uv, tex);

}

vec4 renderMain(){

	if(PASSINDEX == 0){
        return renderMainImage();
	}
    else if (PASSINDEX == 1) {
        
        vec4 img = texture(buffA, _uv);

        if (traceMix < 0.1) {
            return img;
        }

        vec2 uwu = ( ( _uv - 0.5 ) / ( 1.0 + (zoom_direction*mix(1.0, syn_BassLevel, bass_zoom)) ) + 0.5 );
        float size = 0.01;
        float flow = sin(bass_time);
        uwu.x += img.r*size + flow;
        uwu.y += img.g*size + flow;
        uwu -= img.b*size + flow;
        vec4 fragColor = vec4(0.0);

        float thresh = 0.0;       
        if(img.x <= thresh && img.y <= thresh && img.z <= thresh) {
            img = mix(img, texture(syn_FinalPass, uwu), traceMix);
        }

        return img;
        

    }
    if(PASSINDEX == 2.0){

        vec4 img = texture(postFXPass, _uv);
        vec2 uv = _uv;
        vec2 uvL = _uv;

        if (uv.x < 0.5) {
            uv.x = 1.0 - uv.x;
        }

        if (uvL.x > 0.5) {
            uvL.x = 1.0 - uvL.x;
        }

		vec4 mirImg = texture(postFXPass, uv);
        vec4 mirImgL = texture(postFXPass, uvL);

        return mix(mix(img, mirImgL, audio_mir), mirImg, mix(l_vert_mirr_mix, syn_ToggleOnBeat, audio_mir));
	}

    if(PASSINDEX == 3.0){
        
        vec4 img = texture(rVertMirrPass, _uv);
        vec2 uv = _uv;

        if (uv.y < 0.5) {
            uv.y = 1.0 - uv.y;
        }

		vec4 mirImg = texture(rVertMirrPass, uv);

        return mix(img, mirImg, hor_mirr_mix);
	}
}