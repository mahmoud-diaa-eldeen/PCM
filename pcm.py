import sys, io, os
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.append(r"C:\Program Files\Lumerical\v241\api\python")
import lumapi

# ================== الإعدادات ==================
HEAT_FILE = r"E:\optics\Lumerical\PCM\PCM_Heat.ldev"
FDTD_FILE = r"E:\optics\Lumerical\PCM\PCM.fsp"

# لازم الطول الموجي يكون داخل مدى المصدر/المونيتور في PCM.fsp
# المشروع الحالي: source = 400..700 nm، monitor "T" = 400..700 nm
WL       = 550e-9
MAT_A    = "Sb2S3_Amorphous"  # أسماء المواد في FDTD
MAT_C    = "Sb2S3_Crystalline"
PCM_NAME = "PCM_Patch"
MON_NAME = "T"               # اسم المونيتور الفعلي في PCM.fsp (مش T_monitor)
MON_PATH = "::model::" + MON_NAME
MODE_MON = "mode_expansion"  # مونيتور حساب neff (هيتضاف تلقائيًا لو مش موجود)

# نبضة التسخين: 200 ns تشغيل + 10 ns هبوط + تبريد لحد 400 ns
# ده اللي بيحدد الزمن الكلي للـ HEAT. لو عايز تبريد أطول زوّد T_END.
T_ON   = 200e-9
T_FALL = 10e-9
T_END  = 400e-9
PULSE_TABLE = np.array([[0.0], [1e-9], [T_ON], [T_ON + T_FALL], [T_END]])

# قدرة التسخين الفعلية بالوات (مهم: هي اللي في "transient value table" مش خاصية power)
P_HEATER = 0.02

# قيم Xf اللي هنشغّل عندها FDTD (كل نقطة ~80 ثانية). زوّد 0.25/0.75 لو عايز منحنى أدق.
XF_TARGETS = [0.0, 0.5, 1.0]

# neff بيتحسب من مونيتور mode expansion جوه نفس ملف FDTD،
# وده مش محتاج FDTD run كامل. الـ FDTD run مطلوب فقط لـ T_optical.
# لو الذاكرة مش كفاية سيبها False؛ neff هيتحسب عادي.
RUN_FDTD = False

# حدود PCM_Patch (متر)
PCM_BOX = dict(x1=-1021.55e-9, x2=978.452e-9,
               y1=-550e-9,     y2=-50e-9,
               z1=400e-9,      z2=420e-9)

# حركيات JMA (Shafiee: GST)
kB   = 8.617333e-5   # eV/K
Ea   = 2.2           # eV
nu   = 1.79e24       # 1/s
m    = 3.0           # Avrami
Tg   = 453.0         # K
Tm   = 890.0         # K

# ================== 1) HEAT ==================
print("-> تشغيل HEAT...")
heat = lumapi.DEVICE(filename=HEAT_FILE)
heat.switchtolayout()

# مهم جدًا:
# 1) حفظ T فقط يمنع خطأ reshape في lumapi ويقلل حجم الملف جدًا
# 2) جدول transient time steps هو اللي بيحدد زمن النبضة والزمن الكلي
heat.getObjectById("::model::HEAT")["spatial results"] = "T"
heat.getObjectById("::model::HEAT")["transient min time step"] = 5e-9  # يقلل الذاكرة والزمن
power_bc = heat.getObjectById("::model::HEAT::boundary conditions::power")
power_bc["transient time steps"] = PULSE_TABLE
# القيمة الفعلية للقدرة محفوظة في transient value table
# (خاصية power وحدها لا تغيّر الحل في وضع transient)
power_bc["transient value table"] = np.array([[0.0], [P_HEATER], [P_HEATER], [0.0], [0.0]])
power_bc["power"] = P_HEATER

# مسح نتائج أي run قديم عشان نضمن إن الحل الجديد هيتحسب بالظبط
heat.eval("switchtolayout; save;")

heat.run()

try:
    # الصح هو "thermal" مش "temperature"
    thermal = heat.getresult("HEAT", "thermal")

    T_matrix = np.array(thermal["T"])
    if T_matrix.ndim == 3:
        T_matrix = T_matrix[:, :, 0]          # (vertices, time)
    T_matrix = np.squeeze(T_matrix)

    t_s = np.array(thermal["t"]).ravel()      # الزمن بالثواني
    x   = np.array(thermal["x"]).ravel()
    y   = np.array(thermal["y"]).ravel()
    z   = np.array(thermal["z"]).ravel()

    # تطبيق القناع (Mask) لطبقة PCM
    b = PCM_BOX
    mask = (x >= b['x1']) & (x <= b['x2']) & \
           (y >= b['y1']) & (y <= b['y2']) & \
           (z >= b['z1']) & (z <= b['z2'])

    if np.any(mask):
        T_pcm = np.max(T_matrix[mask, :], axis=0)
        T_avg = np.mean(T_matrix[mask, :], axis=0)
        ok = 1
    else:
        print("تنبيه: لم يتم العثور على نقاط داخل حدود PCM_BOX، سيتم استخدام الحرارة العظمى.")
        T_pcm = np.max(T_matrix, axis=0)
        T_avg = np.mean(T_matrix, axis=0)
        ok = 0

except Exception as e:
    print(f"حدث خطأ أثناء استخراج البيانات من HEAT: {e}")
    heat.close()
    sys.exit(1)

heat.close()

print("عدد نقاط الزمن:", len(t_s))
print("t (ns):", np.round(t_s * 1e9, 2))
print("T_pcm max (K):", np.round(T_pcm, 1))
print("T_avg max (K):", round(float(np.max(T_avg)), 1))

# ================== 2) JMA تراكمي ==================
def rate(T):
    T = np.asarray(T, float)
    R = nu * np.exp(-Ea / (kB * T))
    return np.where((T > Tg) & (T < Tm), R, 0.0)

dt  = np.diff(t_s, prepend=0.0)
Xf  = np.zeros_like(t_s)
tau = 0.0
for i in range(len(t_s)):
    if T_pcm[i] >= Tm:
        tau = 0.0
    else:
        tau += rate(T_pcm[i]) * dt[i]
    Xf[i] = 1 - np.exp(-tau**m)
print("Xf max:", float(np.max(Xf)), "| Xf final:", float(Xf[-1]))

# ================== 3) FDTD / neff ==================
# neff من mode expansion monitor جوه نفس ملف FDTD (سريع، من غير run كامل).
# T_optical فقط هو اللي محتاج FDTD run كامل.

# نحاول نقرأ T_optical من CSV قديم لو موجود، عشان لو RUN_FDTD=False
# منضيّعش نتائج الـ FDTD السابقة.
OUT_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pcm_results.csv")
OLD_T = None
if os.path.exists(OUT_CSV):
    try:
        old = np.genfromtxt(OUT_CSV, delimiter=",", names=True)
        if old.size == len(t_s) and old.dtype.names and "T_optical" in old.dtype.names:
            old_t    = np.atleast_1d(old["t_s"]).astype(float)
            old_xf   = np.atleast_1d(old["Xf"]).astype(float)
            old_topt = np.atleast_1d(old["T_optical"]).astype(float)
            same_grid = (len(old_t) == len(t_s)
                         and np.allclose(old_t, t_s, rtol=1e-6, atol=1e-12)
                         and np.allclose(old_xf, Xf, rtol=1e-4, atol=1e-5))
            if (not same_grid) or (not np.any(np.isfinite(old_topt))):
                OLD_T = None
            else:
                OLD_T = old_topt
                print("تم العثور على T_optical من CSV سابق — هنستخدمه لو RUN_FDTD=False.")
    except Exception as e:
        print("تنبيه: تعذر قراءة CSV القديم:", e)

print("\n-> تجهيز FDTD...")
fdtd = lumapi.FDTD(filename=FDTD_FILE)
fdtd.switchtolayout()
f0 = 299792458.0 / WL

def get_nk(mat):
    v = np.array(fdtd.getfdtdindex(mat, f0, 0.99*f0, 1.01*f0)).flatten()[0]
    return v.real, abs(v.imag)

n_a, k_a = get_nk(MAT_A)
n_c, k_c = get_nk(MAT_C)
print(f"amorphous: n={n_a:.4f} k={k_a:.5f} | crystalline: n={n_c:.4f} k={k_c:.5f}")

def mix_nk(x):
    ea = (n_a + 1j*k_a)**2; ec = (n_c + 1j*k_c)**2
    L  = (1-x)*(ea-1)/(ea+2) + x*(ec-1)/(ec+2)
    nk = np.sqrt((1+2*L)/(1-L))
    return nk.real, abs(nk.imag)

def find_mode_monitor(fdtd):
    for cand in ("::model::" + MODE_MON, "::" + MODE_MON):
        try:
            if str(fdtd.getnamed(cand, "type")) == "ModeExpansionMonitor":
                return cand
        except Exception:
            pass
    return None

def ensure_mode_monitor(fdtd, wl):
    """يتأكد إن فيه ModeExpansionMonitor عند x=0 جوه PCM_Patch، ويضبطه على wl."""
    fdtd.switchtolayout()
    obj = find_mode_monitor(fdtd)
    if obj is None:
        # نسجّل الكائنات قبل الإضافة عشان نعرف اسم/مسار المونيتور الجديد
        fdtd.groupscope("::")
        fdtd.selectall()
        before = set(fdtd.getid().split("\n"))
        fdtd.eval("addmodeexpansion;")
        fdtd.groupscope("::")
        fdtd.selectall()
        after = set(fdtd.getid().split("\n"))
        new = [x for x in (after - before) if x]
        if not new:
            raise RuntimeError("addmodeexpansion لم ينشئ مونيتور")
        new_obj = new[0]
        parent = "::".join(new_obj.split("::")[:-1])
        fdtd.setnamed(new_obj, "name", MODE_MON)
        obj = f"{parent}::{MODE_MON}"
        print(f"تم إضافة مونيتور {MODE_MON} تلقائيًا")

    # إعدادات المونيتور: 2D X-normal عند x=0 جوه PCM_Patch
    fdtd.setnamed(obj, "monitor type", "2D X-normal")
    fdtd.setnamed(obj, "x", 0.0)
    fdtd.setnamed(obj, "y", -0.25e-6)
    fdtd.setnamed(obj, "y span", 2.4e-6)
    fdtd.setnamed(obj, "z", 0.0)
    fdtd.setnamed(obj, "z span", 2.3e-6)
    fdtd.setnamed(obj, "mode selection", "fundamental mode")
    fdtd.setnamed(obj, "use source limits", 0)
    fdtd.setnamed(obj, "wavelength center", wl)
    fdtd.setnamed(obj, "wavelength span", 0.0)
    fdtd.setnamed(obj, "frequency points", 1)
    # ملاحظة: "x span" مش active في 2D X-normal، فمينفعش نضبطها.
    return obj

MIX = "PCM_mix"
fdtd.switchtolayout()

# لو المادة موجودة بالفعل (وهي delete-protected) نحدّثها،
# ولو مش موجودة ننشئها مرة واحدة. مينفعش نعمل delete/add كل مرة.
try:
    info = str(fdtd.getmaterial(MIX))
    exists = "not available" not in info
except Exception:
    exists = False

if exists:
    print(f"المادة {MIX} موجودة — سيتم تحديث n/k")
else:
    fdtd.eval(f'newmat = addmaterial("(n,k) Material"); setmaterial(newmat, "name", "{MIX}");')
    print(f"تم إنشاء المادة {MIX}")

fdtd.setnamed(PCM_NAME, "material", MIX)

# مونيتور neff: يتأكد إنه موجود ويضبطه على 550 nm
try:
    MODE_MON_PATH = ensure_mode_monitor(fdtd, WL)
except Exception as e:
    print("خطأ في إعداد مونيتور neff:", e)
    fdtd.close()
    sys.exit(1)

# رجّع الـ group scope الافتراضي بتاع المشروع (المونيتورات تحت ::model)
try:
    fdtd.groupscope("::model")
except Exception:
    pass

try:
    fdtd.save()
except Exception as e:
    print("تنبيه: لم يتم حفظ PCM.fsp:", e)

print(f"\n{'الحالة':<14}{'n':>8}{'k':>10}{'alpha(1/m)':>14}")
for name, x in [("amorphous", 0.0), ("60% cryst", 0.6), ("crystalline", 1.0)]:
    n, k = mix_nk(x)
    print(f"{name:<14}{n:>8.4f}{k:>10.5f}{4*np.pi*k/WL:>14.3e}")

xf_max = float(np.max(Xf))
if xf_max < 1e-3:
    targets = [0.0]
else:
    targets = sorted(set([0.0] + [t for t in XF_TARGETS if 0.0 < t < xf_max] + [xf_max]))

print("\n-> حساب neff عند قيم Xf:", ", ".join(f"{t:.3f}" for t in targets))
neff_targets = []
T_targets = []
warned_wl = False
for xt in targets:
    n, k = mix_nk(xt)
    fdtd.switchtolayout()
    fdtd.setmaterial(MIX, "Refractive Index", n)
    fdtd.setmaterial(MIX, "Imaginary Refractive Index", k)

    # ---- neff: mode solve سريع على مقطع الـ FDTD (مش محتاج run كامل) ----
    try:
        fdtd.eval(f"select('{MODE_MON_PATH}'); updatemodes(1);")
        d_neff = fdtd.getresult(MODE_MON_PATH, "neff")
        neff = complex(np.array(d_neff["neff"]).ravel()[0])
        lam_mode = float(np.array(d_neff["lambda"]).ravel()[0])
    except Exception as e:
        print(f"خطأ في حساب neff عند Xf={xt:.3f}: {e}")
        fdtd.close()
        sys.exit(1)

    if abs(lam_mode - WL) > 1e-9:
        print(f"تحذير: neff اتحسب عند {lam_mode*1e9:.2f} nm مش {WL*1e9:.2f} nm")

    neff_targets.append(neff)

    if RUN_FDTD:
        fdtd.run()
        d = fdtd.getresult(MON_PATH, "T")
        lam = np.array(d["lambda"]).ravel()
        Tvals = np.real(np.array(d["T"])).ravel()
        if not warned_wl and not (lam.min() <= WL <= lam.max()):
            print(f"تحذير: الطول الموجي {WL*1e9:.1f} nm خارج مدى المونيتور "
                  f"({lam.min()*1e9:.0f}-{lam.max()*1e9:.0f} nm).")
            warned_wl = True
        Tt = float(np.interp(WL, lam, Tvals))
        T_targets.append(Tt)
        print(f"  Xf={xt:.3f} -> n={n:.3f} k={k:.4f}  "
              f"neff={neff.real:.5f}{neff.imag:+.5f}j  "
              f"T_opt({WL*1e9:.0f}nm)={Tt:.4f}")
    else:
        T_targets.append(np.nan)
        tnote = "T_opt من CSV قديم" if OLD_T is not None else "T_opt متخطي"
        print(f"  Xf={xt:.3f} -> n={n:.3f} k={k:.4f}  "
              f"neff={neff.real:.5f}{neff.imag:+.5f}j  ({tnote})")

neff_real_t = np.array([z.real for z in neff_targets])
neff_imag_t = np.array([z.imag for z in neff_targets])

if len(targets) == 1:
    res_neff      = np.full_like(Xf, neff_real_t[0], dtype=float)
    res_neff_imag = np.full_like(Xf, neff_imag_t[0], dtype=float)
    res_T         = np.full_like(Xf, T_targets[0], dtype=float)
else:
    res_neff      = np.interp(Xf, targets, neff_real_t)
    res_neff_imag = np.interp(Xf, targets, neff_imag_t)
    res_T         = np.interp(Xf, targets, T_targets)

# لو RUN_FDTD=False استخدم T_optical القديم لو موجود ومتطابق
if (not RUN_FDTD) and (OLD_T is not None):
    res_T = np.array(OLD_T, dtype=float).copy()
    print("تم استخدام T_optical من CSV سابق (RUN_FDTD=False).")

print("\nملخص neff:")
T_show = np.interp(targets, Xf, res_T) if np.any(np.isfinite(res_T)) else np.array(T_targets)
for xt, z, Tt in zip(targets, neff_targets, T_show):
    tstr = f"{Tt:.4f}" if np.isfinite(Tt) else "—"
    print(f"  Xf={xt:.3f}  neff={z.real:.6f}{z.imag:+.6f}j  T_opt={tstr}")

np.savetxt(OUT_CSV, np.c_[t_s, T_pcm, Xf, res_neff, res_neff_imag, res_T],
           delimiter=",",
           header="t_s,T_pcm_K,Xf,neff_real,neff_imag,T_optical",
           comments="")
print("\nتم. النتائج في:", OUT_CSV)

# ================== 4) رسم النتائج ==================
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t_ns = t_s * 1e9
    png = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pcm_results.png")
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))

    ax = axes[0, 0]
    ax.plot(t_ns, T_pcm, color="tab:red")
    ax.set_xlabel("time (ns)"); ax.set_ylabel("T_PCM (K)")
    ax.set_title("Temperature"); ax.grid(alpha=0.3)

    ax = axes[0, 1]
    ax.plot(t_ns, Xf, color="tab:green")
    ax.set_xlabel("time (ns)"); ax.set_ylabel("Xf")
    ax.set_title("Crystallinity"); ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.plot(t_ns, res_neff, color="tab:blue", label="Re(neff)")
    ax.set_xlabel("time (ns)"); ax.set_ylabel("Re(neff)")
    ax.set_title("Effective index"); ax.grid(alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(t_ns, res_neff_imag, color="tab:orange", alpha=0.7, label="Im(neff)")
    ax2.set_ylabel("Im(neff)")
    l1, lab1 = ax.get_legend_handles_labels()
    l2, lab2 = ax2.get_legend_handles_labels()
    ax.legend(l1 + l2, lab1 + lab2, loc="best")

    ax = axes[1, 1]
    ax.plot(targets, neff_real_t, "o-", color="tab:blue", label="Re(neff)")
    ax.set_xlabel("Xf"); ax.set_ylabel("Re(neff)", color="tab:blue")
    ax.tick_params(axis="y", labelcolor="tab:blue")
    ax.grid(alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(targets, neff_imag_t, "s--", color="tab:orange", label="Im(neff)")
    ax2.set_ylabel("Im(neff)", color="tab:orange")
    ax2.tick_params(axis="y", labelcolor="tab:orange")
    ax.set_title("neff vs Xf")

    fig.suptitle("PCM: HEAT -> Xf -> FDTD neff", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(png, dpi=150)
    plt.close(fig)
    print("تم حفظ الرسم في:", png)
except Exception as e:
    print("تنبيه: تعذر إنشاء الرسم:", e)

try:
    fdtd.close()
except Exception:
    pass
