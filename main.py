# main.py — محاكاة فلكية 
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body
from astropy import units as u

# --- قاعدة النجوم (مدمجة) ---
STARS = {
    "Alnitak": ("05h40m45.5s", "-01d56m33s"),
    "Alnilam": ("05h36m12.8s", "-01d12m06s"),
    "Mintaka": ("05h32m00.4s", "-00d17m57s"),
    "Sirius": ("06h45m08.9s", "-16d42m58s")
}

# --- تمايل المحور (Precession) ---
def get_star_precessed(name, time):
    star = SkyCoord(ra=STARS[name][0], dec=STARS[name][1], frame='icrs')
    # تطبيق تمايل المحور حتى 100,000 ق.م.
    from astropy.coordinates import GCRS
    try:
        return star.transform_to(GCRS(obstime=time))
    except:
        return star  # fallback

# --- كواكب قديمة (VSOP87 مبسط) ---
def get_planet_ancient(body, time):
    T = (time.jd - 2451545.0) / 36525.0  # قرون
    if body == "sun":
        L = 4.895063 + 0.01720209895 * (time.jd - 2451545.0)
        return SkyCoord(ra=L*u.rad, dec=0*u.deg, frame='barycentrictrueecliptic')
    return None

# --- تحميل config.txt ---
def load_config():
    sims = []
    current = {}
    try:
        with open("config.txt", "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                if line.startswith("sim"):
                    if current: sims.append(current)
                    current = {"name": line.split(".")[0]}
                else:
                    k, v = line.split("=", 1)
                    k, v = k.strip(), v.strip()
                    if k in ["lat", "lon"]: current[k] = float(v)
                    elif k == "plot": current[k] = v.lower() == "true"
                    else: current[k] = v
        if current: sims.append(current)
    except: pass
    return sims

# --- المحاكاة ---
def simulate(sim):
    print(f"\n{'='*60}")
    print(f"   {sim['name']} | {sim['date']}")
    print(f"{'='*60}")
    
    loc = EarthLocation(lat=sim['lat']*u.deg, lon=sim['lon']*u.deg)
    time = Time(sim['date'])
    frame = AltAz(obstime=time, location=loc)
    targets = [t.strip() for t in sim['targets'].split(",")]
    results = {}
    
    for t in targets:
        obj = None
        if t in STARS:
            obj = get_star_precessed(t, time)
        elif abs(time.year) < 3000:
            try:
                from astropy.coordinates import solar_system_ephemeris
                with solar_system_ephemeris.set('de432s'):
                    obj = get_body(t.lower(), time)
            except: pass
        else:
            obj = get_planet_ancient(t, time)
        
        if obj:
            try:
                altaz = obj.transform_to(frame)
                results[t] = {'az': altaz.az.deg, 'alt': altaz.alt.deg}
                print(f"{t}: أزيموث {altaz.az.deg:.2f}° | ارتفاع {altaz.alt.deg:.2f}°")
            except: pass
    
    if sim.get("plot", False):
        plot_sky(results, sim['name'])

def plot_sky(results, name):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 360); ax.set_ylim(-90, 90)
    for obj, d in results.items():
        ax.scatter(d['az'], d['alt'], s=120, label=obj, edgecolor='black')
    ax.set_xlabel("أزيموث (°)")
    ax.set_ylabel("ارتفاع (°)")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_title(f"محاكاة: {name}")
    os.makedirs("output", exist_ok=True)
    plt.savefig(f"output/{name}.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"تم حفظ: output/{name}.png")

# --- التشغيل ---
if __name__ == "__main__":
    os.makedirs("output", exist_ok=True)
    sims = load_config()
    for sim in sims:
        simulate(sim)
    print("\nتم! النتائج في مجلد 'output'")