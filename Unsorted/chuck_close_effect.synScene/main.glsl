

#define size 12
#define sampleStep 4

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 m =  mod(fragCoord.xy, vec2(size));
	vec2 uv = fragCoord.xy-m;
    vec4[int(pow(float(size/sampleStep), 2.))] c;// = vec4(0.);//texture(syn_UserImage, uv);
    int x = 0;
    for(int i = 0; i<size; i+=sampleStep){
        for(int j = 0; j<size; j+=sampleStep){
            c[x]=texture(syn_UserImage, (uv+vec2(i,j))/RENDERSIZE.xy);
            x++;
        }
    }
    float t = smoothTimeC*0.1;
    float f = length(m-float(size)*.5)*.666*syn_Intensity-t;
    int i = int(mod(f, float(x)));
    int j = int(mod(f+1., float(x)));
    vec4 o =mix( c[i], c[j], mod(f, 1.));
     
	fragColor = o;
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}