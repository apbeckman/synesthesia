

#define size 12
#define sampleStep 4
float dist(vec2 p0, vec2 pf){return sqrt((pf.x-p0.x)*(pf.x-p0.x)+(pf.y-p0.y)*(pf.y-p0.y));}

float d = dist(RENDERSIZE.xy*0.5,_xy.xy)*(_mouse.x/RENDERSIZE.x+0.1)*0.001;

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 m =  mod(fragCoord.xy, vec2(size));
	vec2 uv = fragCoord.xy-m;
    vec4[int(pow(float(size/sampleStep), 2.))] c;// = vec4(0.);//texture(syn_UserImage, uv);
    int x = 0;
    for(int i = 0; i<size; i+=sampleStep){
        for(int j = 0; j<size; j+=sampleStep){
            c[x]=texture(syn_UserImage, (uv+vec2(i,j)*diff)/RENDERSIZE.xy);
            c[x] = _brightness((c[x]), 1.0+flash*(syn_HighLevel*0.125+0.125*syn_MidHighLevel));
            x++;
        }
    }
    float t = smoothTimeC*0.3;
    float f = length(m-float(size)*.5)*.666*(syn_Intensity*0.5+syn_Presence*0.5)-t;
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