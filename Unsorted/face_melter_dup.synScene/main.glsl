//Face Melter - Psybernautics (Alex Tiemann) - 2022
#define FAR 80.

float smin( float a, float b, float k )
{
    float h = max(k-abs(a-b),0.0);
    return min(a, b) - h*h*0.25/k;
}

vec3 re(vec3 p, float d) {return mod(p - d * .5, d) - d * .5;}

vec2 map(vec3 p) {    

    vec3 q = p;
    vec2 scene = vec2(3.0);
    q.z -= 15.0;
    q.y = re(q, floor(instances)*3.0).y;
    scene.x = length(q.yz) - (0.2);
    return scene;
}

vec2 trace(vec3 ro, vec3 rd){
    vec2 t = vec2(0.0);
    vec2 d = vec2(0.0);

    for (int i = 0; i < 48; i++){
        d = map(ro + rd*t.x);
        if(abs(d.x)<.0005 || t.x>FAR) break;        
        t.x += d.x*.75;
    }
    t.y = d.y;
    return t;
}

// Tetrahedral normal, to save a couple of "map" calls. Courtesy of IQ.
vec3 getNormal( in vec3 p ){
    vec2 e = vec2(.0025, -.0025); 
    return normalize(

        e.xyy * map(p + e.xyy).x + 

        e.yyx * map(p + e.yyx).x + 

        e.yxy * map(p + e.yxy).x + 

        e.xxx * map(p + e.xxx).x);

}

vec3 _grad3(vec3 col1, vec3 col2, vec3 col3, float mixVal){

    mixVal *= 2.0;

    float mix1 = clamp(mixVal,0.0,1.0);

    float mix2 = clamp(mixVal-1.0, 0.0, 1.0);

    return mix(mix(col1, col2, mix1), mix(col2, col3, mix2), step(1.0, mixVal));

}

vec3 getObjectColor(vec3 p, float material){

    vec3 col = palette;
    float sdf = map(p-p*_noise(p)*0.1).x;
    vec3 paletteColor = mix(_grad3(vec3(1.0) - col, col, vec3(1.0 - col), material), _palette(sdf, col, col, col, col), noise_tex);
    return paletteColor*10.0;
}

vec3 doColor(in vec3 sp, in vec3 rd, in vec3 sn, in vec3 lp, vec2 t){
    vec3 ld = lp-sp;
    float lDist = max(length(ld), .001);
    ld /= lDist;
    float atten = 1. / (1. + lDist*.2 + lDist*lDist*.1);
    float diff = max(dot(sn, ld), 0.01);
    float spec = pow(max( dot( reflect(ld, sn), -rd ), 0.1), 2.0);
    vec3 objCol = getObjectColor(sp*2.0, t.y);
    if(syn_MediaType > 0.5) {
        vec3 texCol = _loadUserImage().xyz;
        texCol = pow(texCol, vec3(2.0))*10.0;
        objCol = mix(mix(objCol, texCol, show_media), objCol/texCol, divide_media);
    }
    vec3 sceneCol = (objCol*diff*(syn_MediaType > 0.5 ? mix(50.0, 10.0, show_media) : 50.0) + objCol*spec*20.0) * atten;
    return sceneCol;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    int AA = 1;
    vec4 totalC = vec4(0.0);
    for( int m=0; m<AA; m++ )
    for( int nm=0; nm<AA; nm++ )
        {
            vec2 uv = (fragCoord.xy - RENDERSIZE.xy*.5) / RENDERSIZE.y;
            vec3 ro = vec3(0.0, 0.0, -30.0);
            vec3 lp = ro;
            vec3 lk = ro + vec3(0.0,0.0,15.0);
            float fFOV=3.14159/(4.0);
            vec3 forward=normalize(lk-ro);
            vec3 right=normalize(vec3(forward.z,0.,-forward.x));
            vec3 up=cross(forward,right);
            vec3 rd=normalize(forward+fFOV*uv.x*right+fFOV*uv.y*up);

            //Based on Perspective Camera - circuit bending
            vec2 v = (_xy.xy / RENDERSIZE.xy) + RENDERSIZE.y;
            v.x = fract(v.x);
            v.y -= dynamic_time*0.75;
            v = vec2(_fbm(v*(2.0)));
            float pi = 2.0*PI;
            float twpi = 2.0*PI;
            float th =  v.y * pi, ph = v.x * twpi;
            vec3 ssp = vec3( sin(ph) * cos(th), sin(th), cos(ph) * cos(th) );            
            rd = mix(rd, normalize(ssp),mix(0.1, 0.4, syn_Presence));
            rd = normalize(vec3(rd.xy, sqrt(max(rd.z*rd.z - dot(rd.xy, rd.xy)*.085, 0.))));

            vec2 t = trace(ro, rd);
            t.y = dot(v, v);
            vec3 sp = ro + rd*t.x;
            ro += rd*t.x;
            vec3 sn = getNormal(ro);
            vec3 sceneColor = doColor(ro, rd, sn, lp, t);     
            float fogF = smoothstep(0., .99, t.x/FAR);
            sceneColor = mix(sceneColor, vec3(0.0), fogF); 
            totalC += vec4(sqrt(clamp(sceneColor, 0., 1.)), 1);
    }

    return clamp(totalC, 0.0, 1.0);
 } 

vec2 mirrorCoords(vec2 uvIn){

    if (mod(uvIn.x, 2.0) > 1.0){
        uvIn.x = 1.0-uvIn.x;
    }
    if (mod(uvIn.y, 2.0) > 1.0){
        uvIn.y = 1.0-uvIn.y;
    }
    return mod(uvIn, 1.0);
}



vec2 rotateCenter(vec2 uvIn, float amount){
  uvIn.y += (RENDERSIZE.x-RENDERSIZE.y)/RENDERSIZE.x;
  uvIn*=vec2(1.0, RENDERSIZE.y/RENDERSIZE.x);
  _uv2uvc(uvIn);
  uvIn = _rotate(uvIn, amount);
  _uvc2uv(uvIn);
  uvIn/=vec2(1.0, RENDERSIZE.y/RENDERSIZE.x);
  uvIn.y -= (RENDERSIZE.x-RENDERSIZE.y)/RENDERSIZE.x;
  return uvIn;
}

vec4 renderMain(){
	if(PASSINDEX == 0.0){
        vec4 fragColor = renderMainImage();
		return syn_MediaType > 0.5 && show_media > 0.5 ? clamp(fragColor + fragColor*_loadUserImage(), 0.0, 1.0) : fragColor;
	}
    if (PASSINDEX == 1.0) {
        vec4 img = texture(postFXPass, _uv);
        vec4 temp = vec4(0.0);
        if (traceMix < 0.1) {
            return img;
        }
        vec2 uwu = ( ( _uv - 0.5 ) / ( 1.0 + zoom) + 0.5 );
        _uv2uvc(uwu);
        float size = 0.015;
        float bs = syn_BassHits*0.01;
        float ms = syn_BassHits*0.01;
        float hs = syn_BassHits*0.01;
        vec4 imgfb = (texture(tp, _uv + vec2(size, 0.0)) + texture(tp, _uv + vec2(-size, 0.0)) + texture(tp, _uv + vec2(0.0, size)) +texture(tp, _uv + vec2(0.0, -size)) +
        texture(secondPass, _uv + vec2(size, size)) + texture(secondPass, _uv + vec2(-size, -size)) + texture(secondPass, _uv + vec2(-size, size)) +texture(secondPass, _uv + vec2(size, -size))+
        texture(secondPass, _uv + vec2(bs, hs)) + texture(secondPass, _uv + vec2(-ms, -bs)) + texture(secondPass, _uv + vec2(-hs, ms)) +texture(secondPass, _uv + vec2(bs, -hs)));
        float as = size*0.5;
        size *= 0.75;
        vec3 imgfbhsv = _rgb2hsv(imgfb.xyz/12.0);
        imgfb /= mix(12.0, 10.0, clamp(imgfbhsv.z, 0.0, 1.0));
        imgfb.xyz = _hsv2rgb(imgfb.xyz);
        uwu.x += imgfb.b*size + syn_MidPresence*as;
        uwu.y += imgfb.g*size + syn_BassPresence*as;
        uwu -= imgfb.r*size + syn_HighPresence*as;
        vec4 fragColor = vec4(0.0);
        _uvc2uv(uwu);
        float thresh = mix(0.01, 0.01+syn_BassHits*0.5, heavy_blend);       
        if(img.x <= thresh || img.y <= thresh || img.z <= thresh) {
            vec2 o_texCoords = uwu;
            o_texCoords = rotateCenter(o_texCoords, fbRotate*PI);
            vec4 fb = texture(secondPass, o_texCoords);
            fb.xyz *= 1.01;
            fb.xyz = _rgb2hsv(fb.xyz);
            fb.x += mix(0.0, 0.01*_noise(vec3(uwu,imgfbhsv.x)), rainbow);
            fb.y += mix(0.0, 0.0005, fb.x);
            fb.z = clamp(fb.z, 0.0, 1.0);
            fb.xyz = _hsv2rgb(fb.xyz);
            img = mix(img, fb, traceMix);
        }
        return clamp(img, 0.0, 1.0);
    }
    if (PASSINDEX == 2.0) {
        vec4 img = texture(secondPass, _uv);
        vec2 p = _uv;
        p = mix(p, _uvc, _noise((p*20.0 + vec2(0.0, dynamic_time))));
        return texture(secondPass, mix(_uv, p*_uv, mix(0.0, 0.1+0.5*syn_Presence, distort)));
    }
    if(PASSINDEX == 3.0) {
        vec2 p = _uv;
        p.x = abs(1.0 - p.x);//-sin(bass_time - (length(p*0.0625)-bass_time))*0.1;
        _uv2uvc(p);
        float a = atan(p.x, p.y);
        float r = length(p) + 0.15;
        vec2 uv = vec2( 0.15/r + 0.075*dynamic_time, a/3.1415927 );
        vec2 uv2 = vec2( uv.x, atan(p.y,abs(p.x))/3.1415927 );
        uv = mirrorCoords(uv);
        uv.y = 1.0 - uv.y;
        uv2 = mirrorCoords(uv2);
        uv2.y = 1.0 - uv2.y;
        vec3 col = textureGrad( tp, uv, dFdx(uv2), dFdy(uv2) ).xyz;
        vec4 img = texture(tp, _uv);
        col = col*r;
        vec4 fragColor = mix(img, vec4(col, 1.0), tunnel);
        return fragColor;
    }
    
}