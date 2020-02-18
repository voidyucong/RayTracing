#include <iostream>
#include <math.h>
#include <functional>
#include <random>
inline double clamp(double x, double min, double max){ return x < min ? 0 : x > max ? max : x; } 
inline double Random(float min = 0.f, float max = 1.f) {
	static std::uniform_real_distribution<double> distribution(min, max);
    static std::mt19937 generator;
    static std::function<double()> rand_generator = std::bind(distribution, generator);
    return rand_generator();
}
struct Vec3 {
	Vec3() { val[0] = val[1] = val[2] = 0.f; }
	Vec3(float v) { val[0] = val[1] = val[2] = v; }
	Vec3(float v1, float v2, float v3) { val[0] = v1; val[1] = v2; val[2] = v3; }
	inline Vec3& operator+=(const Vec3& r) { val[0] += r.val[0], val[1] += r.val[1], val[2] += r.val[2]; return *this; }
	inline Vec3& operator*=(float v) { val[0] *= v, val[1] *= v, val[2] *= v; return *this; }
	inline Vec3& operator/=(float n) { val[0] /= n, val[1] /= n, val[2] /= n; return *this; }
	inline Vec3 operator-() const { return Vec3(-val[0], -val[1], -val[2]); }
	inline float operator[](int i) const { return val[i]; }
	inline float length() const { return sqrt(val[0] * val[0] + val[1] * val[1] + val[2] * val[2]); }
	float val[3];
};
inline Vec3 operator+(const Vec3& l, const Vec3& r) { return Vec3(l.val[0] + r.val[0], l.val[1] + r.val[1], l.val[2] + r.val[2]); }
inline Vec3 operator+(const Vec3& l, float n) { return Vec3(l.val[0] + n, l.val[1] + n, l.val[2] + n); }
inline Vec3 operator+(float n, const Vec3& r) { return Vec3(r.val[0] + n, r.val[1] + n, r.val[2] + n); }
inline Vec3 operator-(float n, const Vec3& r) { return Vec3(r.val[0] - n, r.val[1] - n, r.val[2] - n); }
inline Vec3 operator-(const Vec3& l, const Vec3& r) { return Vec3(l.val[0] - r.val[0], l.val[1] - r.val[1], l.val[2] - r.val[2]); }
inline Vec3 operator*(const Vec3& l, const Vec3& r) { return Vec3(l.val[0] * r.val[0], l.val[1] * r.val[1], l.val[2] * r.val[2]); }
inline Vec3 operator*(const Vec3& l, float n) { return Vec3(l.val[0] * n, l.val[1] * n, l.val[2] * n); }
inline Vec3 operator*(float n, const Vec3& r) { return Vec3(r.val[0] * n, r.val[1] * n, r.val[2] * n); }
inline Vec3 operator/(const Vec3& l, float n) { return Vec3(l.val[0] / n, l.val[1] / n, l.val[2] / n); }
inline Vec3 normalize(const Vec3& v) { return v / v.length(); }
inline float dot(const Vec3& l, const Vec3& r) { return l.val[0] * r.val[0] + l.val[1] * r.val[1] + l.val[2] * r.val[2]; }
inline Vec3 cross(const Vec3& l, const Vec3& r) { return Vec3(l.val[1] * r.val[2] - l.val[2] * r.val[1], l.val[2] * r.val[0] - l.val[0] * r.val[2], l.val[0] * r.val[1] - l.val[1] * r.val[0]); }
inline Vec3 reflect(const Vec3& v, const Vec3& n) { return v - n * dot(v, n) * 2; }
inline bool refract(const Vec3& v, const Vec3& n, float ni_over_nt, Vec3& refracted) {
    Vec3 uv = normalize(v);
    float dt = dot(uv, n);
    float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
    if (discriminant > 0) { refracted = (uv - n * dt) * ni_over_nt - n * sqrt(discriminant); return true; }
    return false;
}
inline float schlick(float cosine, float ref_idx) {
    float r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 *= r0;
    return r0 + (1 - r0) * pow((1 - cosine), 5);
}
Vec3 RandomInSphere() { return normalize(2.f * Vec3(Random(), Random(), Random()) - Vec3(1.f)) * Random(); }
struct Material;
struct HitRecord {
	float t;
	Vec3 p, normal;
	Material* material;
};
struct Ray {
	Ray() {}
	Ray(Vec3 o, Vec3 r) : origin(o), dir(normalize(r)) {}
	Vec3 origin, dir;
};
Ray GetRay(float x, float y) { return Ray(Vec3(0.f), Vec3(2.f * x - 1.f, 2.f * y - 1.f, -1.f)); }
struct Object { virtual bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& rec) const = 0; };
struct Sphere : public Object {
	Sphere(const Vec3& c, float r, Material* m) : center(c), radius(r), material(m) {}
	virtual bool Hit(const Ray& ray, float t_min, float t_max, HitRecord& rec) const {
		Vec3 oc = ray.origin - center;
	    float a = dot(ray.dir, ray.dir);
	    float b = 2.0 * dot(oc, ray.dir);
	    float c = dot(oc, oc) - radius * radius;
	    float discriminant = b * b - 4 * a * c;
	    if (discriminant > 0) {
	        float temp = (-b - sqrt(discriminant)) / (2.f * a);
	        if (temp < t_max && temp > t_min) { 
	        	rec.t = temp;
	        	rec.p = ray.origin + rec.t * ray.dir;
	        	rec.normal = (rec.p - center) / radius;
	        	rec.material = material;
	        	return true; 
	        }
	        temp = (-b + sqrt(discriminant)) / (2.f * a);
	        if (temp < t_max && temp > t_min) { 
	        	rec.t = temp;
	        	rec.p = ray.origin + rec.t * ray.dir;
	        	rec.normal = (rec.p - center) / radius;
	        	rec.material = material;
	        	return true; 
	        }
	    }
	    return false;
	}
	Vec3 center;
	float radius;
	Material* material;
};
struct Material {
	virtual bool Scatter(const Ray& ray, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const = 0;
	virtual Vec3 Emitted(const Vec3& p) const { return Vec3(0.f); }
};
struct Lambertian : public Material {
	Lambertian(const Vec3& c) : albedo(c) {}
	virtual bool Scatter(const Ray& ray, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
		scattered = Ray(rec.p, rec.normal + RandomInSphere());
		attenuation = albedo;
		return true;
	}
	Vec3 albedo;
};
struct Metal : public Material {
	Metal(const Vec3& c) : albedo(c) {}
	virtual bool Scatter(const Ray& ray, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
		scattered = Ray(rec.p, reflect(ray.dir, rec.normal));
		attenuation = albedo;
		return true;
	}
	Vec3 albedo;
};
struct Dielectric : public Material {
	Dielectric(const Vec3& c, float ri) : albedo(c), refractiveIndice(ri) {}
	virtual bool Scatter(const Ray& ray, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const {
		Vec3 outward_normal, refracted, reflected = reflect(ray.dir, rec.normal);
		float dt, ni_over_nt, reflect_prob, cosine;
		if ((dt = dot(ray.dir, rec.normal)) > 0) {
			outward_normal = -rec.normal;
			ni_over_nt = refractiveIndice;
			cosine = refractiveIndice * dt / ray.dir.length();
		}
		else {
			outward_normal = rec.normal;
			ni_over_nt = 1.0 / refractiveIndice;
			cosine = -dt / ray.dir.length();
		}
		reflect_prob = refract(ray.dir, outward_normal, ni_over_nt, refracted) ? schlick(cosine, refractiveIndice) : 1.f;
		scattered = Ray(rec.p, Random() < reflect_prob ? reflected : refracted);
		attenuation = albedo;
		return true;
	}
	float refractiveIndice;
	Vec3 albedo;
};
struct DiffLight : public Material {
	DiffLight(const Vec3& c) : albedo(c) {}
	virtual bool Scatter(const Ray& ray, const HitRecord& rec, Vec3& attenuation, Ray& scattered) const { return false; }
	virtual Vec3 Emitted(const Vec3& p) const { return albedo; }
	Vec3 albedo;
};
float w = 600, h = 600;
int spp = 600, max_depth = 30;
#define COL(c) clamp(int(255.f * sqrt(c)), 0.f, 255.f)
std::vector<Object*> world = {
	new Sphere(Vec3(0.f, -1e5 - 2, -2.f), 1e5, new Lambertian(Vec3(0.75))),
	new Sphere(Vec3(0.f, 1e5 + 2, -2.f), 1e5, new Lambertian(Vec3(0.75))),
	new Sphere(Vec3(1e5 + 2, -1.f, -2.f), 1e5, new Lambertian(Vec3(0.75, 0.25, 0.25))),
	new Sphere(Vec3(-1e5 - 2, -1.f, -2.f), 1e5, new Lambertian(Vec3(0.25, 0.25, 0.75))),
	new Sphere(Vec3(0.f, -1.f, -1e5 - 4.f), 1e5, new Lambertian(Vec3(0.75))),
	new Sphere(Vec3(0.7f, -1.f, -2.5f), 0.5f, new Lambertian(Vec3(1.f))), 
	new Sphere(Vec3(-1.f, -1.f, -3.f), 0.5f, new Metal(Vec3(1.f))),
	new Sphere(Vec3(0.f, -1.f, -2.f), 0.5f, new Dielectric(Vec3(1.f), 1.5f)),
	new Sphere(Vec3(0.f, 1.2f, -3.f), 0.3f, new DiffLight(Vec3(11.f))),
};
Vec3 color(const Ray& ray, int depth) {
	HitRecord rec;
	bool hitted = false;
	double closest = 1000;
	for (int i = 0; i < world.size(); ++i) if (world[i]->Hit(ray, 0.001, closest, rec)) { closest = rec.t; hitted = true; }
	if (hitted){
		Vec3 attenuation(1.f), emitted = rec.material->Emitted(rec.p);
		Ray scattered;
		if (depth < max_depth && rec.material->Scatter(ray, rec, attenuation, scattered)) { return emitted + attenuation * color(scattered, depth + 1); }
		return emitted;
	}
	return Vec3(0.f);
}
int main() {
	std::cout << "P3\n" << w << " " << h << "\n255\n";
	for (int j = h - 1; j >= 0; --j) {
		for (int i = 0; i < w; ++i) {
			Vec3 col(0.f);
			for (int s = 0; s < spp; ++s)
				col += color(GetRay(float(i + Random()) / w, float(j + Random()) / h), 1);
			col /= spp;
			std::cout << COL(col[0]) << " " << COL(col[1]) << " " << COL(col[2]) << "\n";
		}
	}
	return 0;
}
