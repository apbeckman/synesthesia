

float hash3d(vec3 p3) {
    //thank you mr. hoskins
   	p3  = fract(p3 * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z); 
}

float valueNoise(vec3 v) {
    vec3 vfloor=floor(v);
 	float aaa=hash3d(vfloor);
    float aab=hash3d(vfloor+vec3(0,0,1));
    float aba=hash3d(vfloor+vec3(0,1,0));
    float abb=hash3d(vfloor+vec3(0,1,1));
    float baa=hash3d(vfloor+vec3(1,0,0));
    float bab=hash3d(vfloor+vec3(1,0,1));
    float bba=hash3d(vfloor+vec3(1,1,0));
    float bbb=hash3d(vfloor+vec3(1,1,1));
    
    vec3 vfrac=v-vfloor;
    vfrac=smoothstep(0.,1.,vfrac);
    
    float aa=mix(aaa,aab,vfrac.z);
    float ab=mix(aba,abb,vfrac.z);
    float ba=mix(baa,bab,vfrac.z);
    float bb=mix(bba,bbb,vfrac.z);
    
    return mix(mix(aa,ab,vfrac.y),mix(ba,bb,vfrac.y),vfrac.x);
}

float sphereSDF(vec3 v, vec3 p, float r) {
    return length(v-p)-r;
}

float boxSDF(vec3 v, vec3 p, vec3 b) {
 	vec3 q=abs(v-p)-b;
	return length(max(q,vec3(0.)))+min(0.,max(q.x,max(q.y,q.z)));
}

float SDF(vec3 v) {
    float decay=clamp(valueNoise(30.*v)*valueNoise(80.*v+sin(20.*v.zxy))-0.3,0.,1.);
    float rounding=0.02*(1.-decay);
    //vec3 vmod = v-vec3(0.2,0.2,3);
    vec3 vmod = mod(v,1.);
    float box1 = boxSDF(vmod, vec3(0.5,0.5,0.5), vec3(.12,.12,.12));
    float box2 = boxSDF(vmod, vec3(0.5,0.5,0.5),vec3(1.,.03,.03));
    float box3 = boxSDF(vmod, vec3(0.5,0.5,0.5),vec3(.03,1.,.03));
    float box4 = boxSDF(vmod, vec3(0.5,0.5,0.5),vec3(.03,.03,1.));
    
    
    float morphtime= smoothTimeC*0.1+(v.x-v.y-v.z)*0.02;
    
    float noisemorph = 1.*(floor(morphtime)+smoothstep(0.3,0.7,fract(morphtime)));
    
    return max(min(box1,min(box2,min(box3,box4))),valueNoise(v+vec3(0.5)+noisemorph)-0.5)-rounding;
    		
}

vec3 normal(vec3 v) {
    const float epsilon = 0.0001;
    float sdf1=SDF(v);
 	return normalize(vec3(SDF(v+vec3(epsilon,0,0))-sdf1,
                SDF(v+vec3(0,epsilon,0))-sdf1,
                SDF(v+vec3(0,0,epsilon))-sdf1));
}

mat2 rotate2d(float theta) {
 	return mat2(cos(theta),-sin(theta),sin(theta),cos(theta));   
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // Normalized pixel coordinates (from -1 to 1)
    vec2 uv = (2.*fragCoord-RENDERSIZE.xy)/RENDERSIZE.x;
	vec3 col = vec3(1);
    const int ITERATIONS = 60;
    
    
    vec3 camera = vec3( smoothTime*0.2,0, smoothTime*0.25);
    vec3 pos = camera;
    vec3 dir = normalize(vec3(uv,1.));
    float theta= smoothTimeC*0.1;
    
    dir.yz=rotate2d(sin( smoothTime*0.1))*dir.yz;
    
    dir.xz=rotate2d(theta)*dir.xz;
                  
    vec3 light=normalize(vec3(-1,1,-2));
    
    float fogdist=13.;
    
    const float surf_threshold = 0.002;
    
    for (int i = 0; i < ITERATIONS; i++) {
    	float dist = SDF(pos);
        if (dist < surf_threshold) {
            vec3 light = normalize(vec3(0,0, smoothTimeB)-pos);
            vec3 norm = normal(pos);
            float luminance = 0.5+0.5*dot(norm,light);
            //col = vec3(valueNoise(pos*10.));
            col = vec3(0.2)+0.2*texture(image9,vec2(fract(pos.xy+vec2(0,pos.z)))).xyz;
            col *= luminance;
            col *= 0.4+0.6*(SDF(pos+norm*0.1)-dist)/0.1;
            
           	//col += length(pos-camera)/fogdist;
            
            float fog=float(i)/float(ITERATIONS);
            
            col += fog*fog;
            break;
        }
        pos += dir*dist;
    }
	//col+=vec3(hash3d(vec3(200.*uv,TIME)))*0.05; //ez dither
	col += smoothstep(0.3,2.0,length(uv)); //vignette
    // Output to screen
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}