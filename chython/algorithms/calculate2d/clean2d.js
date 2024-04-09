! function (t, e) {
    "object" == typeof exports && "object" == typeof module ? module.exports = e() : "function" == typeof define && define.amd ? define([], e) : "object" == typeof exports ? exports.$ = e() : t.$ = e()
}(self, (function () {
    return (() => {
        var t = {
            348: t => {
                class e {
                    static clone(t) {
                        let i = Array.isArray(t) ? Array() : {};
                        for (let r in t) {
                            let n = t[r];
                            "function" == typeof n.clone ? i[r] = n.clone() : i[r] = "object" == typeof n ? e.clone(n) : n
                        }
                        return i
                    }
                    static equals(t, e) {
                        if (t.length !== e.length) return !1;
                        let i = t.slice().sort(),
                            r = e.slice().sort();
                        for (var n = 0; n < i.length; n++)
                            if (i[n] !== r[n]) return !1;
                        return !0
                    }
                    static print(t) {
                        if (0 == t.length) return "";
                        let e = "(";
                        for (let i = 0; i < t.length; i++) e += t[i].id ? t[i].id + ", " : t[i] + ", ";
                        return e = e.substring(0, e.length - 2), e + ")"
                    }
                    static each(t, e) {
                        for (let i = 0; i < t.length; i++) e(t[i])
                    }
                    static get(t, e, i) {
                        for (let r = 0; r < t.length; r++)
                            if (t[r][e] == i) return t[r]
                    }
                    static contains(t, e) {
                        if (e.property || e.func) {
                            if (e.func) {
                                for (let i = 0; i < t.length; i++)
                                    if (e.func(t[i])) return !0
                            } else
                                for (let i = 0; i < t.length; i++)
                                    if (t[i][e.property] == e.value) return !0
                        } else
                            for (let i = 0; i < t.length; i++)
                                if (t[i] == e.value) return !0;
                        return !1
                    }
                    static intersection(t, e) {
                        let i = new Array;
                        for (let r = 0; r < t.length; r++)
                            for (let n = 0; n < e.length; n++) t[r] === e[n] && i.push(t[r]);
                        return i
                    }
                    static unique(t) {
                        let e = {};
                        return t.filter((function (t) {
                            return void 0 === e[t] && (e[t] = !0)
                        }))
                    }
                    static count(t, e) {
                        let i = 0;
                        for (let r = 0; r < t.length; r++) t[r] === e && i++;
                        return i
                    }
                    static toggle(t, e) {
                        let i = Array(),
                            r = !1;
                        for (let n = 0; n < t.length; n++) t[n] !== e ? i.push(t[n]) : r = !0;
                        return r || i.push(e), i
                    }
                    static remove(t, e) {
                        let i = Array();
                        for (let r = 0; r < t.length; r++) t[r] !== e && i.push(t[r]);
                        return i
                    }
                    static removeUnique(t, e) {
                        let i = t.indexOf(e);
                        return i > -1 && t.splice(i, 1), t
                    }
                    static removeAll(t, e) {
                        return t.filter((function (t) {
                            return -1 === e.indexOf(t)
                        }))
                    }
                    static merge(t, e) {
                        let i = new Array(t.length + e.length);
                        for (let e = 0; e < t.length; e++) i[e] = t[e];
                        for (let r = 0; r < e.length; r++) i[t.length + r] = e[r];
                        return i
                    }
                    static containsAll(t, e) {
                        let i = 0;
                        for (let r = 0; r < t.length; r++)
                            for (let n = 0; n < e.length; n++) t[r] === e[n] && i++;
                        return i === e.length
                    }
                    static sortByAtomicNumberDesc(t) {
                        let e = t.map((function (t, e) {
                            return {
                                index: e,
                                value: t.atomicNumber.split(".").map(Number)
                            }
                        }));
                        return e.sort((function (t, e) {
                            let i = Math.min(e.value.length, t.value.length),
                                r = 0;
                            for (; r < i && e.value[r] === t.value[r];) r++;
                            return r === i ? e.value.length - t.value.length : e.value[r] - t.value[r]
                        })), e.map((function (e) {
                            return t[e.index]
                        }))
                    }
                    static deepCopy(t) {
                        let i = Array();
                        for (let r = 0; r < t.length; r++) {
                            let n = t[r];
                            i[r] = n instanceof Array ? e.deepCopy(n) : n
                        }
                        return i
                    }
                }
                t.exports = e
            },
            427: (t, e, i) => {
                const r = i(348);
                i(843), i(421);
                class n {
                    constructor(t, e = "-") {
                        this.element = 1 === t.length ? t.toUpperCase() : t, this.drawExplicit = !1, this.ringbonds = Array(), this.rings = Array(), this.bondType = e, this.branchBond = null, this.isBridge = !1, this.isBridgeNode = !1, this.originalRings = Array(), this.bridgedRing = null, this.anchoredRings = Array(), this.bracket = null, this.plane = 0, this.attachedPseudoElements = {}, this.hasAttachedPseudoElements = !1, this.isDrawn = !0, this.isConnectedToRing = !1, this.neighbouringElements = Array(), this.isPartOfAromaticRing = t !== this.element, this.bondCount = 0, this.chirality = "", this.isStereoCenter = !1, this.priority = 0, this.mainChain = !1, this.hydrogenDirection = "down", this.subtreeDepth = 1, this.hasHydrogen = !1, this.class = void 0
                    }
                    addNeighbouringElement(t) {
                        this.neighbouringElements.push(t)
                    }
                    attachPseudoElement(t, e, i = 0, r = 0) {
                        null === i && (i = 0), null === r && (r = 0);
                        let n = i + t + r;
                        this.attachedPseudoElements[n] ? this.attachedPseudoElements[n].count += 1 : this.attachedPseudoElements[n] = {
                            element: t,
                            count: 1,
                            hydrogenCount: i,
                            previousElement: e,
                            charge: r
                        }, this.hasAttachedPseudoElements = !0
                    }
                    getAttachedPseudoElements() {
                        let t = {},
                            e = this;
                        return Object.keys(this.attachedPseudoElements).sort().forEach((function (i) {
                            t[i] = e.attachedPseudoElements[i]
                        })), t
                    }
                    getAttachedPseudoElementsCount() {
                        return Object.keys(this.attachedPseudoElements).length
                    }
                    isHeteroAtom() {
                        return "C" !== this.element && "H" !== this.element
                    }
                    addAnchoredRing(t) {
                        r.contains(this.anchoredRings, {
                            value: t
                        }) || this.anchoredRings.push(t)
                    }
                    getRingbondCount() {
                        return this.ringbonds.length
                    }
                    backupRings() {
                        this.originalRings = Array(this.rings.length);
                        for (let t = 0; t < this.rings.length; t++) this.originalRings[t] = this.rings[t]
                    }
                    restoreRings() {
                        this.rings = Array(this.originalRings.length);
                        for (let t = 0; t < this.originalRings.length; t++) this.rings[t] = this.originalRings[t]
                    }
                    haveCommonRingbond(t, e) {
                        for (let i = 0; i < t.ringbonds.length; i++)
                            for (let r = 0; r < e.ringbonds.length; r++)
                                if (t.ringbonds[i].id == e.ringbonds[r].id) return !0;
                        return !1
                    }
                    neighbouringElementsEqual(t) {
                        if (t.length !== this.neighbouringElements.length) return !1;
                        t.sort(), this.neighbouringElements.sort();
                        for (var e = 0; e < this.neighbouringElements.length; e++)
                            if (t[e] !== this.neighbouringElements[e]) return !1;
                        return !0
                    }
                    getAtomicNumber() {
                        return n.atomicNumbers[this.element]
                    }
                    getMaxBonds() {
                        return n.maxBonds[this.element]
                    }
                    static get maxBonds() {
                        return {
                            H: 1,
                            C: 4,
                            N: 3,
                            O: 2,
                            P: 3,
                            S: 2,
                            B: 3,
                            F: 1,
                            I: 1,
                            Cl: 1,
                            Br: 1
                        }
                    }
                    static get atomicNumbers() {
                        return {
                            H: 1,
                            He: 2,
                            Li: 3,
                            Be: 4,
                            B: 5,
                            b: 5,
                            C: 6,
                            c: 6,
                            N: 7,
                            n: 7,
                            O: 8,
                            o: 8,
                            F: 9,
                            Ne: 10,
                            Na: 11,
                            Mg: 12,
                            Al: 13,
                            Si: 14,
                            P: 15,
                            p: 15,
                            S: 16,
                            s: 16,
                            Cl: 17,
                            Ar: 18,
                            K: 19,
                            Ca: 20,
                            Sc: 21,
                            Ti: 22,
                            V: 23,
                            Cr: 24,
                            Mn: 25,
                            Fe: 26,
                            Co: 27,
                            Ni: 28,
                            Cu: 29,
                            Zn: 30,
                            Ga: 31,
                            Ge: 32,
                            As: 33,
                            Se: 34,
                            Br: 35,
                            Kr: 36,
                            Rb: 37,
                            Sr: 38,
                            Y: 39,
                            Zr: 40,
                            Nb: 41,
                            Mo: 42,
                            Tc: 43,
                            Ru: 44,
                            Rh: 45,
                            Pd: 46,
                            Ag: 47,
                            Cd: 48,
                            In: 49,
                            Sn: 50,
                            Sb: 51,
                            Te: 52,
                            I: 53,
                            Xe: 54,
                            Cs: 55,
                            Ba: 56,
                            La: 57,
                            Ce: 58,
                            Pr: 59,
                            Nd: 60,
                            Pm: 61,
                            Sm: 62,
                            Eu: 63,
                            Gd: 64,
                            Tb: 65,
                            Dy: 66,
                            Ho: 67,
                            Er: 68,
                            Tm: 69,
                            Yb: 70,
                            Lu: 71,
                            Hf: 72,
                            Ta: 73,
                            W: 74,
                            Re: 75,
                            Os: 76,
                            Ir: 77,
                            Pt: 78,
                            Au: 79,
                            Hg: 80,
                            Tl: 81,
                            Pb: 82,
                            Bi: 83,
                            Po: 84,
                            At: 85,
                            Rn: 86,
                            Fr: 87,
                            Ra: 88,
                            Ac: 89,
                            Th: 90,
                            Pa: 91,
                            U: 92,
                            Np: 93,
                            Pu: 94,
                            Am: 95,
                            Cm: 96,
                            Bk: 97,
                            Cf: 98,
                            Es: 99,
                            Fm: 100,
                            Md: 101,
                            No: 102,
                            Lr: 103,
                            Rf: 104,
                            Db: 105,
                            Sg: 106,
                            Bh: 107,
                            Hs: 108,
                            Mt: 109,
                            Ds: 110,
                            Rg: 111,
                            Cn: 112,
                            Uut: 113,
                            Uuq: 114,
                            Uup: 115,
                            Uuh: 116,
                            Uus: 117,
                            Uuo: 118
                        }
                    }
                    static get mass() {
                        return {
                            H: 1,
                            He: 2,
                            Li: 3,
                            Be: 4,
                            B: 5,
                            b: 5,
                            C: 6,
                            c: 6,
                            N: 7,
                            n: 7,
                            O: 8,
                            o: 8,
                            F: 9,
                            Ne: 10,
                            Na: 11,
                            Mg: 12,
                            Al: 13,
                            Si: 14,
                            P: 15,
                            p: 15,
                            S: 16,
                            s: 16,
                            Cl: 17,
                            Ar: 18,
                            K: 19,
                            Ca: 20,
                            Sc: 21,
                            Ti: 22,
                            V: 23,
                            Cr: 24,
                            Mn: 25,
                            Fe: 26,
                            Co: 27,
                            Ni: 28,
                            Cu: 29,
                            Zn: 30,
                            Ga: 31,
                            Ge: 32,
                            As: 33,
                            Se: 34,
                            Br: 35,
                            Kr: 36,
                            Rb: 37,
                            Sr: 38,
                            Y: 39,
                            Zr: 40,
                            Nb: 41,
                            Mo: 42,
                            Tc: 43,
                            Ru: 44,
                            Rh: 45,
                            Pd: 46,
                            Ag: 47,
                            Cd: 48,
                            In: 49,
                            Sn: 50,
                            Sb: 51,
                            Te: 52,
                            I: 53,
                            Xe: 54,
                            Cs: 55,
                            Ba: 56,
                            La: 57,
                            Ce: 58,
                            Pr: 59,
                            Nd: 60,
                            Pm: 61,
                            Sm: 62,
                            Eu: 63,
                            Gd: 64,
                            Tb: 65,
                            Dy: 66,
                            Ho: 67,
                            Er: 68,
                            Tm: 69,
                            Yb: 70,
                            Lu: 71,
                            Hf: 72,
                            Ta: 73,
                            W: 74,
                            Re: 75,
                            Os: 76,
                            Ir: 77,
                            Pt: 78,
                            Au: 79,
                            Hg: 80,
                            Tl: 81,
                            Pb: 82,
                            Bi: 83,
                            Po: 84,
                            At: 85,
                            Rn: 86,
                            Fr: 87,
                            Ra: 88,
                            Ac: 89,
                            Th: 90,
                            Pa: 91,
                            U: 92,
                            Np: 93,
                            Pu: 94,
                            Am: 95,
                            Cm: 96,
                            Bk: 97,
                            Cf: 98,
                            Es: 99,
                            Fm: 100,
                            Md: 101,
                            No: 102,
                            Lr: 103,
                            Rf: 104,
                            Db: 105,
                            Sg: 106,
                            Bh: 107,
                            Hs: 108,
                            Mt: 109,
                            Ds: 110,
                            Rg: 111,
                            Cn: 112,
                            Uut: 113,
                            Uuq: 114,
                            Uup: 115,
                            Uuh: 116,
                            Uus: 117,
                            Uuo: 118
                        }
                    }
                }
                t.exports = n
            },
            841: (t, e, i) => {
                const r = i(474),
                    n = i(614),
                    {
                        getChargeText: s
                    } = (i(929), i(843), i(421), i(537));
                t.exports = class {
                    constructor(t, e, i) {
                        this.canvas = "string" == typeof t || t instanceof String ? document.getElementById(t) : t, this.ctx = this.canvas.getContext("2d"), this.themeManager = e, this.opts = i, this.drawingWidth = 0, this.drawingHeight = 0, this.offsetX = 0, this.offsetY = 0, this.fontLarge = this.opts.fontSizeLarge + "pt Helvetica, Arial, sans-serif", this.fontSmall = this.opts.fontSizeSmall + "pt Helvetica, Arial, sans-serif", this.updateSize(this.opts.width, this.opts.height), this.ctx.font = this.fontLarge, this.hydrogenWidth = this.ctx.measureText("H").width, this.halfHydrogenWidth = this.hydrogenWidth / 2, this.halfBondThickness = this.opts.bondThickness / 2
                    }
                    updateSize(t, e) {
                        this.devicePixelRatio = window.devicePixelRatio || 1, this.backingStoreRatio = this.ctx.webkitBackingStorePixelRatio || this.ctx.mozBackingStorePixelRatio || this.ctx.msBackingStorePixelRatio || this.ctx.oBackingStorePixelRatio || this.ctx.backingStorePixelRatio || 1, this.ratio = this.devicePixelRatio / this.backingStoreRatio, 1 !== this.ratio ? (this.canvas.width = t * this.ratio, this.canvas.height = e * this.ratio, this.canvas.style.width = t + "px", this.canvas.style.height = e + "px", this.ctx.setTransform(this.ratio, 0, 0, this.ratio, 0, 0)) : (this.canvas.width = t * this.ratio, this.canvas.height = e * this.ratio)
                    }
                    setTheme(t) {
                        this.colors = t
                    }
                    scale(t) {
                        let e = -Number.MAX_VALUE,
                            i = -Number.MAX_VALUE,
                            r = Number.MAX_VALUE,
                            n = Number.MAX_VALUE;
                        for (var s = 0; s < t.length; s++) {
                            if (!t[s].value.isDrawn) continue;
                            let o = t[s].position;
                            e < o.x && (e = o.x), i < o.y && (i = o.y), r > o.x && (r = o.x), n > o.y && (n = o.y)
                        }
                        var o = this.opts.padding;
                        e += o, i += o, r -= o, n -= o, this.drawingWidth = e - r, this.drawingHeight = i - n;
                        var h = this.canvas.offsetWidth / this.drawingWidth,
                            a = this.canvas.offsetHeight / this.drawingHeight,
                            l = h < a ? h : a;
                        this.ctx.scale(l, l), this.offsetX = -r, this.offsetY = -n, h < a ? this.offsetY += this.canvas.offsetHeight / (2 * l) - this.drawingHeight / 2 : this.offsetX += this.canvas.offsetWidth / (2 * l) - this.drawingWidth / 2
                    }
                    reset() {
                        this.ctx.setTransform(1, 0, 0, 1, 0, 0)
                    }
                    getColor(t) {
                        return (t = t.toUpperCase()) in this.colors ? this.colors[t] : this.colors.C
                    }
                    drawCircle(t, e, i, n, s = !0, o = !1, h = "") {
                        let a = this.ctx,
                            l = this.offsetX,
                            g = this.offsetY;
                        a.save(), a.lineWidth = 1.5, a.beginPath(), a.arc(t + l, e + g, i, 0, r.twoPI, !0), a.closePath(), o ? (s ? (a.fillStyle = "#f00", a.fill()) : (a.strokeStyle = "#f00", a.stroke()), this.drawDebugText(t, e, h)) : s ? (a.fillStyle = n, a.fill()) : (a.strokeStyle = n, a.stroke()), a.restore()
                    }
                    drawLine(t, e = !1, i = 1) {
                        let r = this.ctx,
                            n = this.offsetX,
                            s = this.offsetY,
                            o = t.clone().shorten(4),
                            h = o.getLeftVector().clone(),
                            a = o.getRightVector().clone();
                        h.x += n, h.y += s, a.x += n, a.y += s, e || (r.save(), r.globalCompositeOperation = "destination-out", r.beginPath(), r.moveTo(h.x, h.y), r.lineTo(a.x, a.y), r.lineCap = "round", r.lineWidth = this.opts.bondThickness + 1.2, r.strokeStyle = this.themeManager.getColor("BACKGROUND"), r.stroke(), r.globalCompositeOperation = "source-over", r.restore()), h = t.getLeftVector().clone(), a = t.getRightVector().clone(), h.x += n, h.y += s, a.x += n, a.y += s, r.save(), r.beginPath(), r.moveTo(h.x, h.y), r.lineTo(a.x, a.y), r.lineCap = "round", r.lineWidth = this.opts.bondThickness;
                        let l = this.ctx.createLinearGradient(h.x, h.y, a.x, a.y);
                        l.addColorStop(.4, this.themeManager.getColor(t.getLeftElement()) || this.themeManager.getColor("C")), l.addColorStop(.6, this.themeManager.getColor(t.getRightElement()) || this.themeManager.getColor("C")), e && (r.setLineDash([1, 1.5]), r.lineWidth = this.opts.bondThickness / 1.5), i < 1 && (r.globalAlpha = i), r.strokeStyle = l, r.stroke(), r.restore()
                    }
                    drawWedge(t, e = 1) {
                        if (isNaN(t.from.x) || isNaN(t.from.y) || isNaN(t.to.x) || isNaN(t.to.y)) return;
                        let i = this.ctx,
                            r = this.offsetX,
                            s = this.offsetY,
                            o = t.clone().shorten(5),
                            h = o.getLeftVector().clone(),
                            a = o.getRightVector().clone();
                        h.x += r, h.y += s, a.x += r, a.y += s, h = t.getLeftVector().clone(), a = t.getRightVector().clone(), h.x += r, h.y += s, a.x += r, a.y += s, i.save();
                        let l = n.normals(h, a);
                        l[0].normalize(), l[1].normalize();
                        let g = h,
                            d = a;
                        t.getRightChiral() && (g = a, d = h);
                        let u = n.add(g, n.multiplyScalar(l[0], this.halfBondThickness)),
                            c = n.add(d, n.multiplyScalar(l[0], 1.5 + this.halfBondThickness)),
                            p = n.add(d, n.multiplyScalar(l[1], 1.5 + this.halfBondThickness)),
                            f = n.add(g, n.multiplyScalar(l[1], this.halfBondThickness));
                        i.beginPath(), i.moveTo(u.x, u.y), i.lineTo(c.x, c.y), i.lineTo(p.x, p.y), i.lineTo(f.x, f.y);
                        let v = this.ctx.createRadialGradient(a.x, a.y, this.opts.bondLength, a.x, a.y, 0);
                        v.addColorStop(.4, this.themeManager.getColor(t.getLeftElement()) || this.themeManager.getColor("C")), v.addColorStop(.6, this.themeManager.getColor(t.getRightElement()) || this.themeManager.getColor("C")), i.fillStyle = v, i.fill(), i.restore()
                    }
                    drawDashedWedge(t) {
                        if (isNaN(t.from.x) || isNaN(t.from.y) || isNaN(t.to.x) || isNaN(t.to.y)) return;
                        let e = this.ctx,
                            i = this.offsetX,
                            r = this.offsetY,
                            s = t.getLeftVector().clone(),
                            o = t.getRightVector().clone();
                        s.x += i, s.y += r, o.x += i, o.y += r, e.save();
                        let h = n.normals(s, o);
                        h[0].normalize(), h[1].normalize();
                        let a, l, g, d, u = t.getRightChiral(),
                            c = t.clone();
                        u ? (a = o, l = s, c.shortenRight(1), g = c.getRightVector().clone(), d = c.getLeftVector().clone()) : (a = s, l = o, c.shortenLeft(1), g = c.getLeftVector().clone(), d = c.getRightVector().clone()), g.x += i, g.y += r, d.x += i, d.y += r;
                        let p = n.subtract(l, a).normalize();
                        e.strokeStyle = this.themeManager.getColor("C"), e.lineCap = "round", e.lineWidth = this.opts.bondThickness, e.beginPath();
                        let f = t.getLength(),
                            v = 1.25 / (f / (3 * this.opts.bondThickness)),
                            m = !1;
                        for (var b = 0; b < 1; b += v) {
                            let i = n.multiplyScalar(p, b * f),
                                r = n.add(a, i),
                                s = 1.5 * b,
                                o = n.multiplyScalar(h[0], s);
                            !m && b > .5 && (e.stroke(), e.beginPath(), e.strokeStyle = this.themeManager.getColor(t.getRightElement()) || this.themeManager.getColor("C"), m = !0), r.subtract(o), e.moveTo(r.x, r.y), r.add(n.multiplyScalar(o, 2)), e.lineTo(r.x, r.y)
                        }
                        e.stroke(), e.restore()
                    }
                    drawDebugText(t, e, i) {
                        let r = this.ctx;
                        r.save(), r.font = "5px Droid Sans, sans-serif", r.textAlign = "start", r.textBaseline = "top", r.fillStyle = "#ff0000", r.fillText(i, t + this.offsetX, e + this.offsetY), r.restore()
                    }
                    drawBall(t, e, i) {
                        let n = this.ctx;
                        n.save(), n.beginPath(), n.arc(t + this.offsetX, e + this.offsetY, this.opts.bondLength / 4.5, 0, r.twoPI, !1), n.fillStyle = this.themeManager.getColor(i), n.fill(), n.restore()
                    }
                    drawPoint(t, e, i) {
                        let n = this.ctx,
                            s = this.offsetX,
                            o = this.offsetY;
                        n.save(), n.globalCompositeOperation = "destination-out", n.beginPath(), n.arc(t + s, e + o, 1.5, 0, r.twoPI, !0), n.closePath(), n.fill(), n.globalCompositeOperation = "source-over", n.beginPath(), n.arc(t + this.offsetX, e + this.offsetY, .75, 0, r.twoPI, !1), n.fillStyle = this.themeManager.getColor(i), n.fill(), n.restore()
                    }
                    drawText(t, e, i, n, o, h, a, l, g, d = {}) {
                        let u = this.ctx,
                            c = this.offsetX,
                            p = this.offsetY;
                        u.save(), u.textAlign = "start", u.textBaseline = "alphabetic";
                        let f = "",
                            v = 0;
                        a && (f = s(a), u.font = this.fontSmall, v = u.measureText(f).width);
                        let m = "0",
                            b = 0;
                        l > 0 && (m = l.toString(), u.font = this.fontSmall, b = u.measureText(m).width), 1 === a && "N" === i && d.hasOwnProperty("0O") && d.hasOwnProperty("0O-1") && (d = {
                            "0O": {
                                element: "O",
                                count: 2,
                                hydrogenCount: 0,
                                previousElement: "C",
                                charge: ""
                            }
                        }, a = 0), u.font = this.fontLarge, u.fillStyle = this.themeManager.getColor("BACKGROUND");
                        let y = u.measureText(i);
                        y.totalWidth = y.width + v, y.height = parseInt(this.fontLarge, 10);
                        let x = y.width > this.opts.fontSizeLarge ? y.width : this.opts.fontSizeLarge;
                        x /= 1.5, u.globalCompositeOperation = "destination-out", u.beginPath(), u.arc(t + c, e + p, x, 0, r.twoPI, !0), u.closePath(), u.fill(), u.globalCompositeOperation = "source-over";
                        let S = -y.width / 2,
                            A = -y.width / 2;
                        u.fillStyle = this.themeManager.getColor(i), u.fillText(i, t + c + S, e + this.opts.halfFontSizeLarge + p), S += y.width, a && (u.font = this.fontSmall, u.fillText(f, t + c + S, e - this.opts.fifthFontSizeSmall + p), S += v), l > 0 && (u.font = this.fontSmall, u.fillText(m, t + c + A - b, e - this.opts.fifthFontSizeSmall + p), A -= b), u.font = this.fontLarge;
                        let C = 0,
                            R = 0;
                        if (1 === n) {
                            let i = t + c,
                                r = e + p + this.opts.halfFontSizeLarge;
                            C = this.hydrogenWidth, A -= C, "left" === o ? i += A : "right" === o || "up" === o && h || "down" === o && h ? i += S : "up" !== o || h ? "down" !== o || h || (r += this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge, i -= this.halfHydrogenWidth) : (r -= this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge, i -= this.halfHydrogenWidth), u.fillText("H", i, r), S += C
                        } else if (n > 1) {
                            let i = t + c,
                                r = e + p + this.opts.halfFontSizeLarge;
                            C = this.hydrogenWidth, u.font = this.fontSmall, R = u.measureText(n).width, A -= C + R, "left" === o ? i += A : "right" === o || "up" === o && h || "down" === o && h ? i += S : "up" !== o || h ? "down" !== o || h || (r += this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge, i -= this.halfHydrogenWidth) : (r -= this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge, i -= this.halfHydrogenWidth), u.font = this.fontLarge, u.fillText("H", i, r), u.font = this.fontSmall, u.fillText(n, i + this.halfHydrogenWidth + R, r + this.opts.fifthFontSizeSmall), S += C + this.halfHydrogenWidth + R
                        }
                        for (let i in d) {
                            if (!d.hasOwnProperty(i)) continue;
                            let r = 0,
                                n = 0,
                                h = d[i].element,
                                a = d[i].count,
                                l = d[i].hydrogenCount,
                                g = d[i].charge;
                            u.font = this.fontLarge, a > 1 && l > 0 && (r = u.measureText("(").width, n = u.measureText(")").width);
                            let f = u.measureText(h).width,
                                v = 0,
                                m = "",
                                b = 0;
                            C = 0, l > 0 && (C = this.hydrogenWidth), u.font = this.fontSmall, a > 1 && (v = u.measureText(a).width), 0 !== g && (m = s(g), b = u.measureText(m).width), R = 0, l > 1 && (R = u.measureText(l).width), u.font = this.fontLarge;
                            let y = t + c,
                                x = e + p + this.opts.halfFontSizeLarge;
                            u.fillStyle = this.themeManager.getColor(h), a > 0 && (A -= v), a > 1 && l > 0 && ("left" === o ? (A -= n, u.fillText(")", y + A, x)) : (u.fillText("(", y + S, x), S += r)), "left" === o ? (A -= f, u.fillText(h, y + A, x)) : (u.fillText(h, y + S, x), S += f), l > 0 && ("left" === o ? (A -= C + R, u.fillText("H", y + A, x), l > 1 && (u.font = this.fontSmall, u.fillText(l, y + A + C, x + this.opts.fifthFontSizeSmall))) : (u.fillText("H", y + S, x), S += C, l > 1 && (u.font = this.fontSmall, u.fillText(l, y + S, x + this.opts.fifthFontSizeSmall), S += R))), u.font = this.fontLarge, a > 1 && l > 0 && ("left" === o ? (A -= r, u.fillText("(", y + A, x)) : (u.fillText(")", y + S, x), S += n)), u.font = this.fontSmall, a > 1 && ("left" === o ? u.fillText(a, y + A + r + n + C + R + f, x + this.opts.fifthFontSizeSmall) : (u.fillText(a, y + S, x + this.opts.fifthFontSizeSmall), S += v)), 0 !== g && ("left" === o ? u.fillText(m, y + A + r + n + C + R + f, e - this.opts.fifthFontSizeSmall + p) : (u.fillText(m, y + S, e - this.opts.fifthFontSizeSmall + p), S += b))
                        }
                        u.restore()
                    }
                    getChargeText(t) {
                        return 1 === t ? "+" : 2 === t ? "2+" : -1 === t ? "-" : -2 === t ? "2-" : ""
                    }
                    drawDebugPoint(t, e, i = "", r = "#f00") {
                        this.drawCircle(t, e, 2, r, !0, !0, i)
                    }
                    drawAromaticityRing(t) {
                        let e = this.ctx,
                            i = r.apothemFromSideLength(this.opts.bondLength, t.getSize());
                        e.save(), e.strokeStyle = this.themeManager.getColor("C"), e.lineWidth = this.opts.bondThickness, e.beginPath(), e.arc(t.center.x + this.offsetX, t.center.y + this.offsetY, i - this.opts.bondSpacing, 0, 2 * Math.PI, !0), e.closePath(), e.stroke(), e.restore()
                    }
                    clear() {
                        this.ctx.clearRect(0, 0, this.canvas.offsetWidth, this.canvas.offsetHeight)
                    }
                }
            },
            237: (t, e, i) => {
                const r = i(474),
                    n = i(348),
                    s = i(614),
                    o = i(929),
                    h = (i(843), i(826)),
                    a = i(427),
                    l = i(421),
                    g = i(333),
                    d = i(841),
                    u = i(707),
                    c = i(473),
                    p = i(654),
                    f = i(207);
                t.exports = class {
                    constructor(t) {
                        this.graph = null, this.doubleBondConfigCount = 0, this.doubleBondConfig = null, this.ringIdCounter = 0, this.ringConnectionIdCounter = 0, this.canvasWrapper = null, this.totalOverlapScore = 0, this.defaultOptions = {
                            width: 500,
                            height: 500,
                            scale: 0,
                            bondThickness: 1,
                            bondLength: 30,
                            shortBondLength: .8,
                            bondSpacing: .17 * 30,
                            atomVisualization: "default",
                            isomeric: !0,
                            debug: !1,
                            terminalCarbons: !1,
                            explicitHydrogens: !0,
                            overlapSensitivity: .42,
                            overlapResolutionIterations: 1,
                            compactDrawing: !0,
                            fontFamily: "Arial, Helvetica, sans-serif",
                            fontSizeLarge: 11,
                            fontSizeSmall: 3,
                            padding: 10,
                            experimentalSSSR: !1,
                            kkThreshold: .1,
                            kkInnerThreshold: .1,
                            kkMaxIteration: 2e4,
                            kkMaxInnerIteration: 50,
                            kkMaxEnergy: 1e9,
                            themes: {
                                dark: {
                                    C: "#fff",
                                    O: "#e74c3c",
                                    N: "#3498db",
                                    F: "#27ae60",
                                    CL: "#16a085",
                                    BR: "#d35400",
                                    I: "#8e44ad",
                                    P: "#d35400",
                                    S: "#f1c40f",
                                    B: "#e67e22",
                                    SI: "#e67e22",
                                    H: "#aaa",
                                    BACKGROUND: "#141414"
                                },
                                light: {
                                    C: "#222",
                                    O: "#e74c3c",
                                    N: "#3498db",
                                    F: "#27ae60",
                                    CL: "#16a085",
                                    BR: "#d35400",
                                    I: "#8e44ad",
                                    P: "#d35400",
                                    S: "#f1c40f",
                                    B: "#e67e22",
                                    SI: "#e67e22",
                                    H: "#666",
                                    BACKGROUND: "#fff"
                                },
                                oldschool: {
                                    C: "#000",
                                    O: "#000",
                                    N: "#000",
                                    F: "#000",
                                    CL: "#000",
                                    BR: "#000",
                                    I: "#000",
                                    P: "#000",
                                    S: "#000",
                                    B: "#000",
                                    SI: "#000",
                                    H: "#000",
                                    BACKGROUND: "#fff"
                                },
                                solarized: {
                                    C: "#586e75",
                                    O: "#dc322f",
                                    N: "#268bd2",
                                    F: "#859900",
                                    CL: "#16a085",
                                    BR: "#cb4b16",
                                    I: "#6c71c4",
                                    P: "#d33682",
                                    S: "#b58900",
                                    B: "#2aa198",
                                    SI: "#2aa198",
                                    H: "#657b83",
                                    BACKGROUND: "#fff"
                                },
                                "solarized-dark": {
                                    C: "#93a1a1",
                                    O: "#dc322f",
                                    N: "#268bd2",
                                    F: "#859900",
                                    CL: "#16a085",
                                    BR: "#cb4b16",
                                    I: "#6c71c4",
                                    P: "#d33682",
                                    S: "#b58900",
                                    B: "#2aa198",
                                    SI: "#2aa198",
                                    H: "#839496",
                                    BACKGROUND: "#fff"
                                },
                                matrix: {
                                    C: "#678c61",
                                    O: "#2fc079",
                                    N: "#4f7e7e",
                                    F: "#90d762",
                                    CL: "#82d967",
                                    BR: "#23755a",
                                    I: "#409931",
                                    P: "#c1ff8a",
                                    S: "#faff00",
                                    B: "#50b45a",
                                    SI: "#409931",
                                    H: "#426644",
                                    BACKGROUND: "#fff"
                                },
                                github: {
                                    C: "#24292f",
                                    O: "#cf222e",
                                    N: "#0969da",
                                    F: "#2da44e",
                                    CL: "#6fdd8b",
                                    BR: "#bc4c00",
                                    I: "#8250df",
                                    P: "#bf3989",
                                    S: "#d4a72c",
                                    B: "#fb8f44",
                                    SI: "#bc4c00",
                                    H: "#57606a",
                                    BACKGROUND: "#fff"
                                },
                                carbon: {
                                    C: "#161616",
                                    O: "#da1e28",
                                    N: "#0f62fe",
                                    F: "#198038",
                                    CL: "#007d79",
                                    BR: "#fa4d56",
                                    I: "#8a3ffc",
                                    P: "#ff832b",
                                    S: "#f1c21b",
                                    B: "#8a3800",
                                    SI: "#e67e22",
                                    H: "#525252",
                                    BACKGROUND: "#fff"
                                },
                                cyberpunk: {
                                    C: "#ea00d9",
                                    O: "#ff3131",
                                    N: "#0abdc6",
                                    F: "#00ff9f",
                                    CL: "#00fe00",
                                    BR: "#fe9f20",
                                    I: "#ff00ff",
                                    P: "#fe7f00",
                                    S: "#fcee0c",
                                    B: "#ff00ff",
                                    SI: "#ffffff",
                                    H: "#913cb1",
                                    BACKGROUND: "#fff"
                                },
                                gruvbox: {
                                    C: "#665c54",
                                    O: "#cc241d",
                                    N: "#458588",
                                    F: "#98971a",
                                    CL: "#79740e",
                                    BR: "#d65d0e",
                                    I: "#b16286",
                                    P: "#af3a03",
                                    S: "#d79921",
                                    B: "#689d6a",
                                    SI: "#427b58",
                                    H: "#7c6f64",
                                    BACKGROUND: "#fbf1c7"
                                },
                                "gruvbox-dark": {
                                    C: "#ebdbb2",
                                    O: "#cc241d",
                                    N: "#458588",
                                    F: "#98971a",
                                    CL: "#b8bb26",
                                    BR: "#d65d0e",
                                    I: "#b16286",
                                    P: "#fe8019",
                                    S: "#d79921",
                                    B: "#8ec07c",
                                    SI: "#83a598",
                                    H: "#bdae93",
                                    BACKGROUND: "#282828"
                                },
                                custom: {
                                    C: "#222",
                                    O: "#e74c3c",
                                    N: "#3498db",
                                    F: "#27ae60",
                                    CL: "#16a085",
                                    BR: "#d35400",
                                    I: "#8e44ad",
                                    P: "#d35400",
                                    S: "#f1c40f",
                                    B: "#e67e22",
                                    SI: "#e67e22",
                                    H: "#666",
                                    BACKGROUND: "#fff"
                                }
                            }
                        }, this.opts = f.extend(!0, this.defaultOptions, t), this.opts.halfBondSpacing = this.opts.bondSpacing / 2, this.opts.bondLengthSq = this.opts.bondLength * this.opts.bondLength, this.opts.halfFontSizeLarge = this.opts.fontSizeLarge / 2, this.opts.quarterFontSizeLarge = this.opts.fontSizeLarge / 4, this.opts.fifthFontSizeSmall = this.opts.fontSizeSmall / 5, this.theme = this.opts.themes.dark
                    }
                    draw(t, e, i = "light", r = !1) {
                        this.initDraw(t, i, r), this.infoOnly || (this.themeManager = new p(this.opts.themes, i), this.canvasWrapper = new d(e, this.themeManager, this.opts)), r || (this.processGraph(), this.canvasWrapper.scale(this.graph.vertices), this.drawEdges(this.opts.debug), this.drawVertices(this.opts.debug), this.canvasWrapper.reset(), this.opts.debug && (console.log(this.graph), console.log(this.rings), console.log(this.ringConnections)))
                    }
                    edgeRingCount(t) {
                        let e = this.graph.edges[t],
                            i = this.graph.vertices[e.sourceId],
                            r = this.graph.vertices[e.targetId];
                        return Math.min(i.value.rings.length, r.value.rings.length)
                    }
                    getBridgedRings() {
                        let t = Array();
                        for (var e = 0; e < this.rings.length; e++) this.rings[e].isBridged && t.push(this.rings[e]);
                        return t
                    }
                    getFusedRings() {
                        let t = Array();
                        for (var e = 0; e < this.rings.length; e++) this.rings[e].isFused && t.push(this.rings[e]);
                        return t
                    }
                    getSpiros() {
                        let t = Array();
                        for (var e = 0; e < this.rings.length; e++) this.rings[e].isSpiro && t.push(this.rings[e]);
                        return t
                    }
                    printRingInfo() {
                        let t = "";
                        for (var e = 0; e < this.rings.length; e++) {
                            const i = this.rings[e];
                            t += i.id + ";", t += i.members.length + ";", t += i.neighbours.length + ";", t += i.isSpiro ? "true;" : "false;", t += i.isFused ? "true;" : "false;", t += i.isBridged ? "true;" : "false;", t += i.rings.length + ";", t += "\n"
                        }
                        return t
                    }
                    rotateDrawing() {
                        let t = 0,
                            e = 0,
                            i = 0;
                        for (var r = 0; r < this.graph.vertices.length; r++) {
                            let s = this.graph.vertices[r];
                            if (s.value.isDrawn)
                                for (var n = r + 1; n < this.graph.vertices.length; n++) {
                                    let o = this.graph.vertices[n];
                                    if (!o.value.isDrawn) continue;
                                    let h = s.position.distanceSq(o.position);
                                    h > i && (i = h, t = r, e = n)
                                }
                        }
                        let o = -s.subtract(this.graph.vertices[t].position, this.graph.vertices[e].position).angle();
                        if (!isNaN(o)) {
                            let t = o % .523599;
                            for (t < .2617995 ? o -= t : o += .523599 - t, r = 0; r < this.graph.vertices.length; r++) r !== e && this.graph.vertices[r].position.rotateAround(o, this.graph.vertices[e].position);
                            for (r = 0; r < this.rings.length; r++) this.rings[r].center.rotateAround(o, this.graph.vertices[e].position)
                        }
                    }
                    getTotalOverlapScore() {
                        return this.totalOverlapScore
                    }
                    getRingCount() {
                        return this.rings.length
                    }
                    hasBridgedRing() {
                        return this.bridgedRing
                    }
                    getHeavyAtomCount() {
                        let t = 0;
                        for (var e = 0; e < this.graph.vertices.length; e++) "H" !== this.graph.vertices[e].value.element && t++;
                        return t
                    }
                    getMolecularFormula(t = null) {
                        let e = "",
                            i = new Map,
                            r = null === t ? this.graph : new u(t, this.opts.isomeric);
                        for (var n = 0; n < r.vertices.length; n++) {
                            let t = r.vertices[n].value;
                            if (i.has(t.element) ? i.set(t.element, i.get(t.element) + 1) : i.set(t.element, 1), t.bracket && !t.bracket.chirality && (i.has("H") ? i.set("H", i.get("H") + t.bracket.hcount) : i.set("H", t.bracket.hcount)), !t.bracket) {
                                let e = a.maxBonds[t.element] - t.bondCount;
                                t.isPartOfAromaticRing && e--, i.has("H") ? i.set("H", i.get("H") + e) : i.set("H", e)
                            }
                        }
                        if (i.has("C")) {
                            let t = i.get("C");
                            e += "C" + (t > 1 ? t : ""), i.delete("C")
                        }
                        if (i.has("H")) {
                            let t = i.get("H");
                            e += "H" + (t > 1 ? t : ""), i.delete("H")
                        }
                        return Object.keys(a.atomicNumbers).sort().map((t => {
                            if (i.has(t)) {
                                let r = i.get(t);
                                e += t + (r > 1 ? r : "")
                            }
                        })), e
                    }
                    getRingbondType(t, e) {
                        if (t.value.getRingbondCount() < 1 || e.value.getRingbondCount() < 1) return null;
                        for (var i = 0; i < t.value.ringbonds.length; i++)
                            for (var r = 0; r < e.value.ringbonds.length; r++)
                                if (t.value.ringbonds[i].id === e.value.ringbonds[r].id) return "-" === t.value.ringbonds[i].bondType ? e.value.ringbonds[r].bond : t.value.ringbonds[i].bond;
                        return null
                    }
                    initDraw(t, e, i, r) {
                        this.data = t, this.infoOnly = i, this.ringIdCounter = 0, this.ringConnectionIdCounter = 0, this.graph = new u(t, this.opts.isomeric), this.rings = Array(), this.ringConnections = Array(), this.originalRings = Array(), this.originalRingConnections = Array(), this.bridgedRing = !1, this.doubleBondConfigCount = null, this.doubleBondConfig = null, this.highlight_atoms = r, this.initRings(), this.initHydrogens()
                    }
                    processGraph() {
                        this.position(), this.restoreRingInformation(), this.resolvePrimaryOverlaps();
                        let t = this.getOverlapScore();
                        this.totalOverlapScore = this.getOverlapScore().total;
                        for (var e = 0; e < this.opts.overlapResolutionIterations; e++)
                            for (var i = 0; i < this.graph.edges.length; i++) {
                                let e = this.graph.edges[i];
                                if (this.isEdgeRotatable(e)) {
                                    let i = this.graph.getTreeDepth(e.sourceId, e.targetId),
                                        n = this.graph.getTreeDepth(e.targetId, e.sourceId),
                                        s = e.targetId,
                                        o = e.sourceId;
                                    if (i > n && (s = e.sourceId, o = e.targetId), this.getSubtreeOverlapScore(o, s, t.vertexScores).value > this.opts.overlapSensitivity) {
                                        let e = this.graph.vertices[s],
                                            i = this.graph.vertices[o],
                                            n = i.getNeighbours(s);
                                        if (1 === n.length) {
                                            let t = this.graph.vertices[n[0]],
                                                s = t.position.getRotateAwayFromAngle(e.position, i.position, r.toRad(120));
                                            this.rotateSubtree(t.id, i.id, s, i.position);
                                            let o = this.getOverlapScore().total;
                                            o > this.totalOverlapScore ? this.rotateSubtree(t.id, i.id, -s, i.position) : this.totalOverlapScore = o
                                        } else if (2 === n.length) {
                                            if (0 !== i.value.rings.length && 0 !== e.value.rings.length) continue;
                                            let t = this.graph.vertices[n[0]],
                                                s = this.graph.vertices[n[1]];
                                            if (1 === t.value.rings.length && 1 === s.value.rings.length) {
                                                if (t.value.rings[0] !== s.value.rings[0]) continue
                                            } else {
                                                if (0 !== t.value.rings.length || 0 !== s.value.rings.length) continue; {
                                                    let n = t.position.getRotateAwayFromAngle(e.position, i.position, r.toRad(120)),
                                                        o = s.position.getRotateAwayFromAngle(e.position, i.position, r.toRad(120));
                                                    this.rotateSubtree(t.id, i.id, n, i.position), this.rotateSubtree(s.id, i.id, o, i.position);
                                                    let h = this.getOverlapScore().total;
                                                    h > this.totalOverlapScore ? (this.rotateSubtree(t.id, i.id, -n, i.position), this.rotateSubtree(s.id, i.id, -o, i.position)) : this.totalOverlapScore = h
                                                }
                                            }
                                        }
                                        t = this.getOverlapScore()
                                    }
                                }
                            }
                        this.resolveSecondaryOverlaps(t.scores), this.opts.isomeric && this.annotateStereochemistry(), this.opts.compactDrawing && "default" === this.opts.atomVisualization && this.initPseudoElements(), this.rotateDrawing()
                    }
                    initRings() {
                        let t = new Map;
                        for (var e = this.graph.vertices.length - 1; e >= 0; e--) {
                            let r = this.graph.vertices[e];
                            if (0 !== r.value.ringbonds.length)
                                for (var i = 0; i < r.value.ringbonds.length; i++) {
                                    let e = r.value.ringbonds[i].id,
                                        n = r.value.ringbonds[i].bond;
                                    if (t.has(e)) {
                                        let s = r.id,
                                            o = t.get(e)[0],
                                            a = t.get(e)[1],
                                            l = new h(s, o, 1);
                                        l.setBondType(a || n || "-");
                                        let g = this.graph.addEdge(l),
                                            d = this.graph.vertices[o];
                                        r.addRingbondChild(o, i), r.value.addNeighbouringElement(d.value.element), d.addRingbondChild(s, i), d.value.addNeighbouringElement(r.value.element), r.edges.push(g), d.edges.push(g), t.delete(e)
                                    } else t.set(e, [r.id, n])
                                }
                        }
                        let r = c.getRings(this.graph, this.opts.experimentalSSSR);
                        if (null !== r) {
                            for (e = 0; e < r.length; e++) {
                                let t = [...r[e]],
                                    n = this.addRing(new l(t));
                                for (i = 0; i < t.length; i++) this.graph.vertices[t[i]].value.rings.push(n)
                            }
                            for (e = 0; e < this.rings.length - 1; e++)
                                for (i = e + 1; i < this.rings.length; i++) {
                                    let t = this.rings[e],
                                        r = this.rings[i],
                                        n = new g(t, r);
                                    n.vertices.size > 0 && this.addRingConnection(n)
                                }
                            for (e = 0; e < this.rings.length; e++) {
                                let t = this.rings[e];
                                t.neighbours = g.getNeighbours(this.ringConnections, t.id)
                            }
                            for (e = 0; e < this.rings.length; e++) {
                                let t = this.rings[e];
                                this.graph.vertices[t.members[0]].value.addAnchoredRing(t.id)
                            }
                            for (this.backupRingInformation(); this.rings.length > 0;) {
                                let t = -1;
                                for (e = 0; e < this.rings.length; e++) {
                                    let i = this.rings[e];
                                    this.isPartOfBridgedRing(i.id) && !i.isBridged && (t = i.id)
                                }
                                if (-1 === t) break;
                                let i = this.getRing(t),
                                    r = this.getBridgedRingRings(i.id);
                                for (this.bridgedRing = !0, this.createBridgedRing(r, i.members[0]), e = 0; e < r.length; e++) this.removeRing(r[e])
                            }
                        }
                    }
                    initHydrogens() {
                        if (!this.opts.explicitHydrogens)
                            for (var t = 0; t < this.graph.vertices.length; t++) {
                                let e = this.graph.vertices[t];
                                if ("H" !== e.value.element) continue;
                                let i = this.graph.vertices[e.neighbours[0]];
                                i.value.hasHydrogen = !0, (!i.value.isStereoCenter || i.value.rings.length < 2 && !i.value.bridgedRing || i.value.bridgedRing && i.value.originalRings.length < 2) && (e.value.isDrawn = !1)
                            }
                    }
                    getBridgedRingRings(t) {
                        let e = Array(),
                            i = this,
                            r = function (t) {
                                let n = i.getRing(t);
                                e.push(t);
                                for (var s = 0; s < n.neighbours.length; s++) {
                                    let o = n.neighbours[s]; - 1 === e.indexOf(o) && o !== t && g.isBridge(i.ringConnections, i.graph.vertices, t, o) && r(o)
                                }
                            };
                        return r(t), n.unique(e)
                    }
                    isPartOfBridgedRing(t) {
                        for (var e = 0; e < this.ringConnections.length; e++)
                            if (this.ringConnections[e].containsRing(t) && this.ringConnections[e].isBridge(this.graph.vertices)) return !0;
                        return !1
                    }
                    createBridgedRing(t, e) {
                        let i = new Set,
                            r = new Set,
                            s = new Set;
                        for (var o = 0; o < t.length; o++) {
                            let e = this.getRing(t[o]);
                            e.isPartOfBridged = !0;
                            for (var h = 0; h < e.members.length; h++) r.add(e.members[h]);
                            for (h = 0; h < e.neighbours.length; h++) {
                                let i = e.neighbours[h]; - 1 === t.indexOf(i) && s.add(e.neighbours[h])
                            }
                        }
                        let a = new Set;
                        for (let e of r) {
                            let r = this.graph.vertices[e],
                                s = n.intersection(t, r.value.rings);
                            1 === r.value.rings.length || 1 === s.length ? i.add(r.id) : a.add(r.id)
                        }
                        Array();
                        let g = Array();
                        for (let t of a) {
                            let e = this.graph.vertices[t],
                                r = !1;
                            for (let t = 0; t < e.edges.length; t++) 1 === this.edgeRingCount(e.edges[t]) && (r = !0);
                            r ? (e.value.isBridgeNode = !0, i.add(e.id)) : (e.value.isBridge = !0, i.add(e.id))
                        }
                        let d = new l([...i]);
                        for (this.addRing(d), d.isBridged = !0, d.neighbours = [...s], o = 0; o < t.length; o++) d.rings.push(this.getRing(t[o]).clone());
                        for (o = 0; o < d.members.length; o++) this.graph.vertices[d.members[o]].value.bridgedRing = d.id;
                        for (o = 0; o < g.length; o++) this.graph.vertices[g[o]].value.rings = Array();
                        for (let e of i) {
                            let i = this.graph.vertices[e];
                            i.value.rings = n.removeAll(i.value.rings, t), i.value.rings.push(d.id)
                        }
                        for (o = 0; o < t.length; o++)
                            for (h = o + 1; h < t.length; h++) this.removeRingConnectionsBetween(t[o], t[h]);
                        for (let e of s) {
                            let i = this.getRingConnections(e, t);
                            for (h = 0; h < i.length; h++) this.getRingConnection(i[h]).updateOther(d.id, e);
                            this.getRing(e).neighbours.push(d.id)
                        }
                        return d
                    }
                    areVerticesInSameRing(t, e) {
                        for (var i = 0; i < t.value.rings.length; i++)
                            for (var r = 0; r < e.value.rings.length; r++)
                                if (t.value.rings[i] === e.value.rings[r]) return !0;
                        return !1
                    }
                    getCommonRings(t, e) {
                        let i = Array();
                        for (var r = 0; r < t.value.rings.length; r++)
                            for (var n = 0; n < e.value.rings.length; n++) t.value.rings[r] == e.value.rings[n] && i.push(t.value.rings[r]);
                        return i
                    }
                    getLargestOrAromaticCommonRing(t, e) {
                        let i = this.getCommonRings(t, e),
                            r = 0,
                            n = null;
                        for (var s = 0; s < i.length; s++) {
                            let t = this.getRing(i[s]),
                                e = t.getSize();
                            if (t.isBenzeneLike(this.graph.vertices)) return t;
                            e > r && (r = e, n = t)
                        }
                        return n
                    }
                    getVerticesAt(t, e, i) {
                        let r = Array();
                        for (var n = 0; n < this.graph.vertices.length; n++) {
                            let s = this.graph.vertices[n];
                            s.id !== i && s.positioned && t.distanceSq(s.position) <= e * e && r.push(s.id)
                        }
                        return r
                    }
                    getClosestVertex(t) {
                        let e = 99999,
                            i = null;
                        for (var r = 0; r < this.graph.vertices.length; r++) {
                            let n = this.graph.vertices[r];
                            if (n.id === t.id) continue;
                            let s = t.position.distanceSq(n.position);
                            s < e && (e = s, i = n)
                        }
                        return i
                    }
                    addRing(t) {
                        return t.id = this.ringIdCounter++, this.rings.push(t), t.id
                    }
                    removeRing(t) {
                        this.rings = this.rings.filter((function (e) {
                            return e.id !== t
                        })), this.ringConnections = this.ringConnections.filter((function (e) {
                            return e.firstRingId !== t && e.secondRingId !== t
                        }));
                        for (var e = 0; e < this.rings.length; e++) {
                            let i = this.rings[e];
                            i.neighbours = i.neighbours.filter((function (e) {
                                return e !== t
                            }))
                        }
                    }
                    getRing(t) {
                        for (var e = 0; e < this.rings.length; e++)
                            if (this.rings[e].id == t) return this.rings[e]
                    }
                    addRingConnection(t) {
                        return t.id = this.ringConnectionIdCounter++, this.ringConnections.push(t), t.id
                    }
                    removeRingConnection(t) {
                        this.ringConnections = this.ringConnections.filter((function (e) {
                            return e.id !== t
                        }))
                    }
                    removeRingConnectionsBetween(t, e) {
                        let i = Array();
                        for (var r = 0; r < this.ringConnections.length; r++) {
                            let n = this.ringConnections[r];
                            (n.firstRingId === t && n.secondRingId === e || n.firstRingId === e && n.secondRingId === t) && i.push(n.id)
                        }
                        for (r = 0; r < i.length; r++) this.removeRingConnection(i[r])
                    }
                    getRingConnection(t) {
                        for (var e = 0; e < this.ringConnections.length; e++)
                            if (this.ringConnections[e].id == t) return this.ringConnections[e]
                    }
                    getRingConnections(t, e) {
                        let i = Array();
                        for (var r = 0; r < this.ringConnections.length; r++) {
                            let s = this.ringConnections[r];
                            for (var n = 0; n < e.length; n++) {
                                let r = e[n];
                                (s.firstRingId === t && s.secondRingId === r || s.firstRingId === r && s.secondRingId === t) && i.push(s.id)
                            }
                        }
                        return i
                    }
                    getOverlapScore() {
                        let t = 0,
                            e = new Float32Array(this.graph.vertices.length);
                        for (var i = 0; i < this.graph.vertices.length; i++) e[i] = 0;
                        for (i = 0; i < this.graph.vertices.length; i++)
                            for (var r = this.graph.vertices.length; --r > i;) {
                                let n = this.graph.vertices[i],
                                    o = this.graph.vertices[r];
                                if (!n.value.isDrawn || !o.value.isDrawn) continue;
                                let h = s.subtract(n.position, o.position).lengthSq();
                                if (h < this.opts.bondLengthSq) {
                                    let n = (this.opts.bondLength - Math.sqrt(h)) / this.opts.bondLength;
                                    t += n, e[i] += n, e[r] += n
                                }
                            }
                        let n = Array();
                        for (i = 0; i < this.graph.vertices.length; i++) n.push({
                            id: i,
                            score: e[i]
                        });
                        return n.sort((function (t, e) {
                            return e.score - t.score
                        })), {
                            total: t,
                            scores: n,
                            vertexScores: e
                        }
                    }
                    chooseSide(t, e, i) {
                        let r = t.getNeighbours(e.id),
                            s = e.getNeighbours(t.id),
                            o = r.length,
                            h = s.length,
                            a = n.merge(r, s),
                            l = [0, 0];
                        for (var g = 0; g < a.length; g++) this.graph.vertices[a[g]].position.sameSideAs(t.position, e.position, i[0]) ? l[0]++ : l[1]++;
                        let d = [0, 0];
                        for (g = 0; g < this.graph.vertices.length; g++) this.graph.vertices[g].position.sameSideAs(t.position, e.position, i[0]) ? d[0]++ : d[1]++;
                        return {
                            totalSideCount: d,
                            totalPosition: d[0] > d[1] ? 0 : 1,
                            sideCount: l,
                            position: l[0] > l[1] ? 0 : 1,
                            anCount: o,
                            bnCount: h
                        }
                    }
                    setRingCenter(t) {
                        let e = t.getSize(),
                            i = new s(0, 0);
                        for (var r = 0; r < e; r++) i.add(this.graph.vertices[t.members[r]].position);
                        t.center = i.divide(e)
                    }
                    getSubringCenter(t, e) {
                        let i = e.value.originalRings,
                            r = t.center,
                            n = Number.MAX_VALUE;
                        for (var s = 0; s < i.length; s++)
                            for (var o = 0; o < t.rings.length; o++) i[s] === t.rings[o].id && t.rings[o].getSize() < n && (r = t.rings[o].center, n = t.rings[o].getSize());
                        return r
                    }
                    drawEdges(t) {
                        let e = this,
                            i = Array(this.graph.edges.length);
                        if (i.fill(!1), this.graph.traverseBF(0, (function (r) {
                            let n = e.graph.getEdges(r.id);
                            for (var s = 0; s < n.length; s++) {
                                let r = n[s];
                                i[r] || (i[r] = !0, e.drawEdge(r, t))
                            }
                        })), !this.bridgedRing)
                            for (var r = 0; r < this.rings.length; r++) {
                                let t = this.rings[r];
                                this.isRingAromatic(t) && this.canvasWrapper.drawAromaticityRing(t)
                            }
                    }
                    drawEdge(t, e) {
                        let i = this,
                            r = this.graph.edges[t],
                            h = this.graph.vertices[r.sourceId],
                            a = this.graph.vertices[r.targetId],
                            l = h.value.element,
                            g = a.value.element;
                        if (!(h.value.isDrawn && a.value.isDrawn || "default" !== this.opts.atomVisualization)) return;
                        let d = h.position,
                            u = a.position,
                            c = this.getEdgeNormals(r),
                            p = n.clone(c);
                        if (p[0].multiplyScalar(10).add(d), p[1].multiplyScalar(10).add(d), "=" === r.bondType || "=" === this.getRingbondType(h, a) || r.isPartOfAromaticRing && this.bridgedRing) {
                            let t = this.areVerticesInSameRing(h, a),
                                e = this.chooseSide(h, a, p);
                            if (t) {
                                let t = this.getLargestOrAromaticCommonRing(h, a).center;
                                c[0].multiplyScalar(i.opts.bondSpacing), c[1].multiplyScalar(i.opts.bondSpacing);
                                let e = null;
                                e = t.sameSideAs(h.position, a.position, s.add(d, c[0])) ? new o(s.add(d, c[0]), s.add(u, c[0]), l, g) : new o(s.add(d, c[1]), s.add(u, c[1]), l, g), e.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength), r.isPartOfAromaticRing ? this.canvasWrapper.drawLine(e, !0) : this.canvasWrapper.drawLine(e), this.canvasWrapper.drawLine(new o(d, u, l, g))
                            } else if (r.center || h.isTerminal() && a.isTerminal()) {
                                c[0].multiplyScalar(i.opts.halfBondSpacing), c[1].multiplyScalar(i.opts.halfBondSpacing);
                                let t = new o(s.add(d, c[0]), s.add(u, c[0]), l, g),
                                    e = new o(s.add(d, c[1]), s.add(u, c[1]), l, g);
                                this.canvasWrapper.drawLine(t), this.canvasWrapper.drawLine(e)
                            } else if (0 == e.anCount && e.bnCount > 1 || 0 == e.bnCount && e.anCount > 1) {
                                c[0].multiplyScalar(i.opts.halfBondSpacing), c[1].multiplyScalar(i.opts.halfBondSpacing);
                                let t = new o(s.add(d, c[0]), s.add(u, c[0]), l, g),
                                    e = new o(s.add(d, c[1]), s.add(u, c[1]), l, g);
                                this.canvasWrapper.drawLine(t), this.canvasWrapper.drawLine(e)
                            } else if (e.sideCount[0] > e.sideCount[1]) {
                                c[0].multiplyScalar(i.opts.bondSpacing), c[1].multiplyScalar(i.opts.bondSpacing);
                                let t = new o(s.add(d, c[0]), s.add(u, c[0]), l, g);
                                t.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength), this.canvasWrapper.drawLine(t), this.canvasWrapper.drawLine(new o(d, u, l, g))
                            } else if (e.sideCount[0] < e.sideCount[1]) {
                                c[0].multiplyScalar(i.opts.bondSpacing), c[1].multiplyScalar(i.opts.bondSpacing);
                                let t = new o(s.add(d, c[1]), s.add(u, c[1]), l, g);
                                t.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength), this.canvasWrapper.drawLine(t), this.canvasWrapper.drawLine(new o(d, u, l, g))
                            } else if (e.totalSideCount[0] > e.totalSideCount[1]) {
                                c[0].multiplyScalar(i.opts.bondSpacing), c[1].multiplyScalar(i.opts.bondSpacing);
                                let t = new o(s.add(d, c[0]), s.add(u, c[0]), l, g);
                                t.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength), this.canvasWrapper.drawLine(t), this.canvasWrapper.drawLine(new o(d, u, l, g))
                            } else if (e.totalSideCount[0] <= e.totalSideCount[1]) {
                                c[0].multiplyScalar(i.opts.bondSpacing), c[1].multiplyScalar(i.opts.bondSpacing);
                                let t = new o(s.add(d, c[1]), s.add(u, c[1]), l, g);
                                t.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength), this.canvasWrapper.drawLine(t), this.canvasWrapper.drawLine(new o(d, u, l, g))
                            }
                        } else if ("#" === r.bondType) {
                            c[0].multiplyScalar(i.opts.bondSpacing / 1.5), c[1].multiplyScalar(i.opts.bondSpacing / 1.5);
                            let t = new o(s.add(d, c[0]), s.add(u, c[0]), l, g),
                                e = new o(s.add(d, c[1]), s.add(u, c[1]), l, g);
                            this.canvasWrapper.drawLine(t), this.canvasWrapper.drawLine(e), this.canvasWrapper.drawLine(new o(d, u, l, g))
                        } else if ("." === r.bondType);
                        else {
                            let t = h.value.isStereoCenter,
                                e = a.value.isStereoCenter;
                            "up" === r.wedge ? this.canvasWrapper.drawWedge(new o(d, u, l, g, t, e)) : "down" === r.wedge ? this.canvasWrapper.drawDashedWedge(new o(d, u, l, g, t, e)) : this.canvasWrapper.drawLine(new o(d, u, l, g, t, e))
                        }
                        if (e) {
                            let e = s.midpoint(d, u);
                            this.canvasWrapper.drawDebugText(e.x, e.y, "e: " + t)
                        }
                    }
                    drawVertices(t) {
                        var e = this.graph.vertices.length;
                        for (e = 0; e < this.graph.vertices.length; e++) {
                            let i = this.graph.vertices[e],
                                r = i.value,
                                o = 0,
                                h = 0,
                                l = i.value.bondCount,
                                g = r.element,
                                d = a.maxBonds[g] - l,
                                u = i.getTextDirection(this.graph.vertices),
                                c = !(!this.opts.terminalCarbons && "C" === g && !r.hasAttachedPseudoElements) && i.isTerminal(),
                                p = "C" === r.element;
                            if ("N" === r.element && r.isPartOfAromaticRing && (d = 0), r.bracket && (d = r.bracket.hcount, o = r.bracket.charge, h = r.bracket.isotope), "allballs" === this.opts.atomVisualization) this.canvasWrapper.drawBall(i.position.x, i.position.y, g);
                            else if (r.isDrawn && (!p || r.drawExplicit || c || r.hasAttachedPseudoElements) || 1 === this.graph.vertices.length) "default" === this.opts.atomVisualization ? this.canvasWrapper.drawText(i.position.x, i.position.y, g, d, u, c, o, h, this.graph.vertices.length, r.getAttachedPseudoElements()) : "balls" === this.opts.atomVisualization && this.canvasWrapper.drawBall(i.position.x, i.position.y, g);
                            else if (2 === i.getNeighbourCount() && 1 == i.forcePositioned) {
                                let t = this.graph.vertices[i.neighbours[0]].position,
                                    e = this.graph.vertices[i.neighbours[1]].position,
                                    r = s.threePointangle(i.position, t, e);
                                Math.abs(Math.PI - r) < .1 && this.canvasWrapper.drawPoint(i.position.x, i.position.y, g)
                            }
                            if (t) {
                                let t = "v: " + i.id + " " + n.print(r.ringbonds);
                                this.canvasWrapper.drawDebugText(i.position.x, i.position.y, t)
                            }
                        }
                        if (this.opts.debug)
                            for (e = 0; e < this.rings.length; e++) {
                                let t = this.rings[e].center;
                                this.canvasWrapper.drawDebugPoint(t.x, t.y, "r: " + this.rings[e].id)
                            }
                    }
                    position() {
                        let t = null;
                        for (var e = 0; e < this.graph.vertices.length; e++)
                            if (null !== this.graph.vertices[e].value.bridgedRing) {
                                t = this.graph.vertices[e];
                                break
                            } if (null === t)
                            for (e = 0; e < this.rings.length; e++)
                                if (this.rings[e].isBridged) {
                                    t = this.graph.vertices[this.rings[e].members[0]];
                                    break
                                } null === t && this.rings.length > 0 && (t = this.graph.vertices[this.rings[0].members[0]]), null === t && (t = this.graph.vertices[0]), this.createNextBond(t, null, 0)
                    }
                    backupRingInformation() {
                        this.originalRings = Array(), this.originalRingConnections = Array();
                        for (var t = 0; t < this.rings.length; t++) this.originalRings.push(this.rings[t]);
                        for (t = 0; t < this.ringConnections.length; t++) this.originalRingConnections.push(this.ringConnections[t]);
                        for (t = 0; t < this.graph.vertices.length; t++) this.graph.vertices[t].value.backupRings()
                    }
                    restoreRingInformation() {
                        let t = this.getBridgedRings();
                        this.rings = Array(), this.ringConnections = Array();
                        for (var e = 0; e < t.length; e++) {
                            let r = t[e];
                            for (var i = 0; i < r.rings.length; i++) {
                                let t = r.rings[i];
                                this.originalRings[t.id].center = t.center
                            }
                        }
                        for (e = 0; e < this.originalRings.length; e++) this.rings.push(this.originalRings[e]);
                        for (e = 0; e < this.originalRingConnections.length; e++) this.ringConnections.push(this.originalRingConnections[e]);
                        for (e = 0; e < this.graph.vertices.length; e++) this.graph.vertices[e].value.restoreRings()
                    }
                    createRing(t, e = null, i = null, n = null) {
                        if (t.positioned) return;
                        e = e || new s(0, 0);
                        let o = t.getOrderedNeighbours(this.ringConnections),
                            h = i ? s.subtract(i.position, e).angle() : 0,
                            a = r.polyCircumradius(this.opts.bondLength, t.getSize()),
                            l = r.centralAngle(t.getSize());
                        t.centralAngle = l;
                        let d = h,
                            u = this,
                            c = i ? i.id : null;
                        if (-1 === t.members.indexOf(c) && (i && (i.positioned = !1), c = t.members[0]), t.isBridged) {
                            this.graph.kkLayout(t.members.slice(), e, i.id, t, this.opts.bondLength, this.opts.kkThreshold, this.opts.kkInnerThreshold, this.opts.kkMaxIteration, this.opts.kkMaxInnerIteration, this.opts.kkMaxEnergy), t.positioned = !0, this.setRingCenter(t), e = t.center;
                            for (var p = 0; p < t.rings.length; p++) this.setRingCenter(t.rings[p])
                        } else t.eachMember(this.graph.vertices, (function (i) {
                            let r = u.graph.vertices[i];
                            r.positioned || r.setPosition(e.x + Math.cos(d) * a, e.y + Math.sin(d) * a), d += l, (!t.isBridged || t.rings.length < 3) && (r.angle = d, r.positioned = !0)
                        }), c, n ? n.id : null);
                        for (t.positioned = !0, t.center = e, p = 0; p < o.length; p++) {
                            let i = this.getRing(o[p].neighbour);
                            if (i.positioned) continue;
                            let n = g.getVertices(this.ringConnections, t.id, i.id);
                            if (2 === n.length) {
                                t.isFused = !0, i.isFused = !0;
                                let o = this.graph.vertices[n[0]],
                                    h = this.graph.vertices[n[1]],
                                    a = s.midpoint(o.position, h.position),
                                    l = s.normals(o.position, h.position);
                                l[0].normalize(), l[1].normalize();
                                let g = r.polyCircumradius(this.opts.bondLength, i.getSize()),
                                    d = r.apothem(g, i.getSize());
                                l[0].multiplyScalar(d).add(a), l[1].multiplyScalar(d).add(a);
                                let u = l[0];
                                s.subtract(e, l[1]).lengthSq() > s.subtract(e, l[0]).lengthSq() && (u = l[1]);
                                let c = s.subtract(o.position, u),
                                    p = s.subtract(h.position, u); - 1 === c.clockwise(p) ? i.positioned || this.createRing(i, u, o, h) : i.positioned || this.createRing(i, u, h, o)
                            } else if (1 === n.length) {
                                t.isSpiro = !0, i.isSpiro = !0;
                                let o = this.graph.vertices[n[0]],
                                    h = s.subtract(e, o.position);
                                h.invert(), h.normalize();
                                let a = r.polyCircumradius(this.opts.bondLength, i.getSize());
                                h.multiplyScalar(a), h.add(o.position), i.positioned || this.createRing(i, h, o)
                            }
                        }
                        for (p = 0; p < t.members.length; p++) {
                            let e = this.graph.vertices[t.members[p]],
                                i = e.neighbours;
                            for (var f = 0; f < i.length; f++) {
                                let t = this.graph.vertices[i[f]];
                                t.positioned || (t.value.isConnectedToRing = !0, this.createNextBond(t, e, 0))
                            }
                        }
                    }
                    rotateSubtree(t, e, i, r) {
                        let n = this;
                        this.graph.traverseTree(t, e, (function (t) {
                            t.position.rotateAround(i, r);
                            for (var e = 0; e < t.value.anchoredRings.length; e++) {
                                let s = n.rings[t.value.anchoredRings[e]];
                                s && s.center.rotateAround(i, r)
                            }
                        }))
                    }
                    getSubtreeOverlapScore(t, e, i) {
                        let r = this,
                            n = 0,
                            o = new s(0, 0),
                            h = 0;
                        return this.graph.traverseTree(t, e, (function (t) {
                            if (!t.value.isDrawn) return;
                            let e = i[t.id];
                            e > r.opts.overlapSensitivity && (n += e, h++);
                            let s = r.graph.vertices[t.id].position.clone();
                            s.multiplyScalar(e), o.add(s)
                        })), o.divide(n), {
                            value: n / h,
                            center: o
                        }
                    }
                    getCurrentCenterOfMass() {
                        let t = new s(0, 0),
                            e = 0;
                        for (var i = 0; i < this.graph.vertices.length; i++) {
                            let r = this.graph.vertices[i];
                            r.positioned && (t.add(r.position), e++)
                        }
                        return t.divide(e)
                    }
                    getCurrentCenterOfMassInNeigbourhood(t, e = 2 * this.opts.bondLength) {
                        let i = new s(0, 0),
                            r = 0,
                            n = e * e;
                        for (var o = 0; o < this.graph.vertices.length; o++) {
                            let e = this.graph.vertices[o];
                            e.positioned && t.distanceSq(e.position) < n && (i.add(e.position), r++)
                        }
                        return i.divide(r)
                    }
                    resolvePrimaryOverlaps() {
                        let t = Array(),
                            e = Array(this.graph.vertices.length);
                        for (var i = 0; i < this.rings.length; i++) {
                            let s = this.rings[i];
                            for (var r = 0; r < s.members.length; r++) {
                                let i = this.graph.vertices[s.members[r]];
                                if (e[i.id]) continue;
                                e[i.id] = !0;
                                let o = this.getNonRingNeighbours(i.id);
                                if (o.length > 1) {
                                    let e = Array();
                                    for (var n = 0; n < i.value.rings.length; n++) e.push(i.value.rings[n]);
                                    t.push({
                                        common: i,
                                        rings: e,
                                        vertices: o
                                    })
                                } else if (1 === o.length && 2 === i.value.rings.length) {
                                    let e = Array();
                                    for (n = 0; n < i.value.rings.length; n++) e.push(i.value.rings[n]);
                                    t.push({
                                        common: i,
                                        rings: e,
                                        vertices: o
                                    })
                                }
                            }
                        }
                        for (i = 0; i < t.length; i++) {
                            let e = t[i];
                            if (2 === e.vertices.length) {
                                let t = e.vertices[0],
                                    i = e.vertices[1];
                                if (!t.value.isDrawn || !i.value.isDrawn) continue;
                                let r = (2 * Math.PI - this.getRing(e.rings[0]).getAngle()) / 6;
                                this.rotateSubtree(t.id, e.common.id, r, e.common.position), this.rotateSubtree(i.id, e.common.id, -r, e.common.position);
                                let n = this.getOverlapScore(),
                                    s = this.getSubtreeOverlapScore(t.id, e.common.id, n.vertexScores),
                                    o = this.getSubtreeOverlapScore(i.id, e.common.id, n.vertexScores),
                                    h = s.value + o.value;
                                this.rotateSubtree(t.id, e.common.id, -2 * r, e.common.position), this.rotateSubtree(i.id, e.common.id, 2 * r, e.common.position), n = this.getOverlapScore(), s = this.getSubtreeOverlapScore(t.id, e.common.id, n.vertexScores), o = this.getSubtreeOverlapScore(i.id, e.common.id, n.vertexScores), s.value + o.value > h && (this.rotateSubtree(t.id, e.common.id, 2 * r, e.common.position), this.rotateSubtree(i.id, e.common.id, -2 * r, e.common.position))
                            } else 1 === e.vertices.length && e.rings.length
                        }
                    }
                    resolveSecondaryOverlaps(t) {
                        for (var e = 0; e < t.length; e++)
                            if (t[e].score > this.opts.overlapSensitivity) {
                                let i = this.graph.vertices[t[e].id];
                                if (i.isTerminal()) {
                                    let t = this.getClosestVertex(i);
                                    if (t) {
                                        let e = null;
                                        e = t.isTerminal() ? 0 === t.id ? this.graph.vertices[1].position : t.previousPosition : 0 === t.id ? this.graph.vertices[1].position : t.position;
                                        let n = 0 === i.id ? this.graph.vertices[1].position : i.previousPosition;
                                        i.position.rotateAwayFrom(e, n, r.toRad(20))
                                    }
                                }
                            }
                    }
                    getLastVertexWithAngle(t) {
                        let e = 0,
                            i = null;
                        for (; !e && t;) i = this.graph.vertices[t], e = i.angle, t = i.parentVertexId;
                        return i
                    }
                    createNextBond(t, e = null, i = 0, o = !1, h = !1) {
                        if (t.positioned && !h) return;
                        let a = !1;
                        if (e) {
                            let i = this.graph.getEdge(t.id, e.id);
                            "/" !== i.bondType && "\\" !== i.bondType || ++this.doubleBondConfigCount % 2 != 1 || null === this.doubleBondConfig && (this.doubleBondConfig = i.bondType, a = !0, null === e.parentVertexId && t.value.branchBond && ("/" === this.doubleBondConfig ? this.doubleBondConfig = "\\" : "\\" === this.doubleBondConfig && (this.doubleBondConfig = "/")))
                        }
                        if (!h)
                            if (e)
                                if (e.value.rings.length > 0) {
                                    let i = e.neighbours,
                                        r = null,
                                        o = new s(0, 0);
                                    if (null === e.value.bridgedRing && e.value.rings.length > 1)
                                        for (var l = 0; l < i.length; l++) {
                                            let t = this.graph.vertices[i[l]];
                                            if (n.containsAll(t.value.rings, e.value.rings)) {
                                                r = t;
                                                break
                                            }
                                        }
                                    if (null === r) {
                                        for (l = 0; l < i.length; l++) {
                                            let t = this.graph.vertices[i[l]];
                                            t.positioned && this.areVerticesInSameRing(t, e) && o.add(s.subtract(t.position, e.position))
                                        }
                                        o.invert().normalize().multiplyScalar(this.opts.bondLength).add(e.position)
                                    } else o = r.position.clone().rotateAround(Math.PI, e.position);
                                    t.previousPosition = e.position, t.setPositionFromVector(o), t.positioned = !0
                                } else {
                                    let r = new s(this.opts.bondLength, 0);
                                    r.rotate(i), r.add(e.position), t.setPositionFromVector(r), t.previousPosition = e.position, t.positioned = !0
                                }
                            else {
                                let e = new s(this.opts.bondLength, 0);
                                e.rotate(r.toRad(-60)), t.previousPosition = e, t.setPosition(this.opts.bondLength, 0), t.angle = r.toRad(-60), null === t.value.bridgedRing && (t.positioned = !0)
                            }
                        if (null !== t.value.bridgedRing) {
                            let e = this.getRing(t.value.bridgedRing);
                            if (!e.positioned) {
                                let i = s.subtract(t.previousPosition, t.position);
                                i.invert(), i.normalize();
                                let n = r.polyCircumradius(this.opts.bondLength, e.members.length);
                                i.multiplyScalar(n), i.add(t.position), this.createRing(e, i, t)
                            }
                        } else if (t.value.rings.length > 0) {
                            let e = this.getRing(t.value.rings[0]);
                            if (!e.positioned) {
                                let i = s.subtract(t.previousPosition, t.position);
                                i.invert(), i.normalize();
                                let n = r.polyCircumradius(this.opts.bondLength, e.getSize());
                                i.multiplyScalar(n), i.add(t.position), this.createRing(e, i, t)
                            }
                        } else {
                            t.value.isStereoCenter;
                            let i = t.getNeighbours(),
                                h = Array();
                            for (l = 0; l < i.length; l++) this.graph.vertices[i[l]].value.isDrawn && h.push(i[l]);
                            e && (h = n.remove(h, e.id));
                            let g = t.getAngle();
                            if (1 === h.length) {
                                let i = this.graph.vertices[h[0]];
                                if ("#" === t.value.bondType || e && "#" === e.value.bondType || "=" === t.value.bondType && e && 0 === e.value.rings.length && "=" === e.value.bondType && "-" !== t.value.branchBond) t.value.drawExplicit = !1, e && (this.graph.getEdge(t.id, e.id).center = !0), this.graph.getEdge(t.id, i.id).center = !0, ("#" === t.value.bondType || e && "#" === e.value.bondType) && (i.angle = 0), i.drawExplicit = !0, this.createNextBond(i, t, g + i.angle);
                                else if (e && e.value.rings.length > 0) {
                                    let e = r.toRad(60),
                                        n = -e,
                                        o = new s(this.opts.bondLength, 0),
                                        h = new s(this.opts.bondLength, 0);
                                    o.rotate(e).add(t.position), h.rotate(n).add(t.position);
                                    let a = this.getCurrentCenterOfMass(),
                                        l = o.distanceSq(a),
                                        d = h.distanceSq(a);
                                    i.angle = l < d ? n : e, this.createNextBond(i, t, g + i.angle)
                                } else {
                                    let r = t.angle;
                                    if (e && e.neighbours.length > 3 ? r = r > 0 ? Math.min(1.0472, r) : r < 0 ? Math.max(-1.0472, r) : 1.0472 : r || (r = this.getLastVertexWithAngle(t.id).angle, r || (r = 1.0472)), e && !a) {
                                        let e = this.graph.getEdge(t.id, i.id).bondType;
                                        "/" === e ? ("/" === this.doubleBondConfig || "\\" === this.doubleBondConfig && (r = -r), this.doubleBondConfig = null) : "\\" === e && ("/" === this.doubleBondConfig ? r = -r : this.doubleBondConfig, this.doubleBondConfig = null)
                                    }
                                    i.angle = o ? r : -r, this.createNextBond(i, t, g + i.angle)
                                }
                            } else if (2 === h.length) {
                                let i = t.angle;
                                i || (i = 1.0472);
                                let r = this.graph.getTreeDepth(h[0], t.id),
                                    n = this.graph.getTreeDepth(h[1], t.id),
                                    s = this.graph.vertices[h[0]],
                                    o = this.graph.vertices[h[1]];
                                s.value.subtreeDepth = r, o.value.subtreeDepth = n;
                                let a = this.graph.getTreeDepth(e ? e.id : null, t.id);
                                e && (e.value.subtreeDepth = a);
                                let l = 0,
                                    d = 1;
                                "C" === o.value.element && "C" !== s.value.element && n > 1 && r < 5 ? (l = 1, d = 0) : "C" !== o.value.element && "C" === s.value.element && r > 1 && n < 5 ? (l = 0, d = 1) : n > r && (l = 1, d = 0);
                                let u = this.graph.vertices[h[l]],
                                    c = this.graph.vertices[h[d]],
                                    p = (this.graph.getEdge(t.id, u.id), this.graph.getEdge(t.id, c.id), !1);
                                a < r && a < n && (p = !0), c.angle = i, u.angle = -i, "\\" === this.doubleBondConfig ? "\\" === c.value.branchBond && (c.angle = -i, u.angle = i) : "/" === this.doubleBondConfig && "/" === c.value.branchBond && (c.angle = -i, u.angle = i), this.createNextBond(c, t, g + c.angle, p), this.createNextBond(u, t, g + u.angle, p)
                            } else if (3 === h.length) {
                                let i = this.graph.getTreeDepth(h[0], t.id),
                                    n = this.graph.getTreeDepth(h[1], t.id),
                                    s = this.graph.getTreeDepth(h[2], t.id),
                                    o = this.graph.vertices[h[0]],
                                    a = this.graph.vertices[h[1]],
                                    l = this.graph.vertices[h[2]];
                                o.value.subtreeDepth = i, a.value.subtreeDepth = n, l.value.subtreeDepth = s, n > i && n > s ? (o = this.graph.vertices[h[1]], a = this.graph.vertices[h[0]], l = this.graph.vertices[h[2]]) : s > i && s > n && (o = this.graph.vertices[h[2]], a = this.graph.vertices[h[0]], l = this.graph.vertices[h[1]]), e && e.value.rings.length < 1 && o.value.rings.length < 1 && a.value.rings.length < 1 && l.value.rings.length < 1 && 1 === this.graph.getTreeDepth(a.id, t.id) && 1 === this.graph.getTreeDepth(l.id, t.id) && this.graph.getTreeDepth(o.id, t.id) > 1 ? (o.angle = -t.angle, t.angle >= 0 ? (a.angle = r.toRad(30), l.angle = r.toRad(90)) : (a.angle = -r.toRad(30), l.angle = -r.toRad(90)), this.createNextBond(o, t, g + o.angle), this.createNextBond(a, t, g + a.angle), this.createNextBond(l, t, g + l.angle)) : (o.angle = 0, a.angle = r.toRad(90), l.angle = -r.toRad(90), this.createNextBond(o, t, g + o.angle), this.createNextBond(a, t, g + a.angle), this.createNextBond(l, t, g + l.angle))
                            } else if (4 === h.length) {
                                let e = this.graph.getTreeDepth(h[0], t.id),
                                    i = this.graph.getTreeDepth(h[1], t.id),
                                    n = this.graph.getTreeDepth(h[2], t.id),
                                    s = this.graph.getTreeDepth(h[3], t.id),
                                    o = this.graph.vertices[h[0]],
                                    a = this.graph.vertices[h[1]],
                                    l = this.graph.vertices[h[2]],
                                    d = this.graph.vertices[h[3]];
                                o.value.subtreeDepth = e, a.value.subtreeDepth = i, l.value.subtreeDepth = n, d.value.subtreeDepth = s, i > e && i > n && i > s ? (o = this.graph.vertices[h[1]], a = this.graph.vertices[h[0]], l = this.graph.vertices[h[2]], d = this.graph.vertices[h[3]]) : n > e && n > i && n > s ? (o = this.graph.vertices[h[2]], a = this.graph.vertices[h[0]], l = this.graph.vertices[h[1]], d = this.graph.vertices[h[3]]) : s > e && s > i && s > n && (o = this.graph.vertices[h[3]], a = this.graph.vertices[h[0]], l = this.graph.vertices[h[1]], d = this.graph.vertices[h[2]]), o.angle = -r.toRad(36), a.angle = r.toRad(36), l.angle = -r.toRad(108), d.angle = r.toRad(108), this.createNextBond(o, t, g + o.angle), this.createNextBond(a, t, g + a.angle), this.createNextBond(l, t, g + l.angle), this.createNextBond(d, t, g + d.angle)
                            }
                        }
                    }
                    getCommonRingbondNeighbour(t) {
                        let e = t.neighbours;
                        for (var i = 0; i < e.length; i++) {
                            let r = this.graph.vertices[e[i]];
                            if (n.containsAll(r.value.rings, t.value.rings)) return r
                        }
                        return null
                    }
                    isPointInRing(t) {
                        for (var e = 0; e < this.rings.length; e++) {
                            let i = this.rings[e];
                            if (!i.positioned) continue;
                            let n = r.polyCircumradius(this.opts.bondLength, i.getSize()),
                                s = n * n;
                            if (t.distanceSq(i.center) < s) return !0
                        }
                        return !1
                    }
                    isEdgeInRing(t) {
                        let e = this.graph.vertices[t.sourceId],
                            i = this.graph.vertices[t.targetId];
                        return this.areVerticesInSameRing(e, i)
                    }
                    isEdgeRotatable(t) {
                        let e = this.graph.vertices[t.sourceId],
                            i = this.graph.vertices[t.targetId];
                        return !("-" !== t.bondType || e.isTerminal() || i.isTerminal() || e.value.rings.length > 0 && i.value.rings.length > 0 && this.areVerticesInSameRing(e, i))
                    }
                    isRingAromatic(t) {
                        for (var e = 0; e < t.members.length; e++)
                            if (!this.graph.vertices[t.members[e]].value.isPartOfAromaticRing) return !1;
                        return !0
                    }
                    getEdgeNormals(t) {
                        let e = this.graph.vertices[t.sourceId].position,
                            i = this.graph.vertices[t.targetId].position;
                        return s.units(e, i)
                    }
                    getNonRingNeighbours(t) {
                        let e = Array(),
                            i = this.graph.vertices[t],
                            r = i.neighbours;
                        for (var s = 0; s < r.length; s++) {
                            let t = this.graph.vertices[r[s]];
                            0 === n.intersection(i.value.rings, t.value.rings).length && 0 == t.value.isBridge && e.push(t)
                        }
                        return e
                    }
                    annotateStereochemistry() {
                        for (var t = 0; t < this.graph.vertices.length; t++) {
                            let s = this.graph.vertices[t];
                            if (!s.value.isStereoCenter) continue;
                            let o = s.getNeighbours(),
                                h = o.length,
                                a = Array(h);
                            for (var e = 0; e < h; e++) {
                                let t = new Uint8Array(this.graph.vertices.length),
                                    r = Array(Array());
                                t[s.id] = 1, this.visitStereochemistry(o[e], s.id, t, r, 10, 0);
                                for (var i = 0; i < r.length; i++) r[i].sort((function (t, e) {
                                    return e - t
                                }));
                                a[e] = [e, r]
                            }
                            let l = 0,
                                g = 0;
                            for (e = 0; e < a.length; e++)
                                for (a[e][1].length > l && (l = a[e][1].length), i = 0; i < a[e][1].length; i++) a[e][1][i].length > g && (g = a[e][1][i].length);
                            for (e = 0; e < a.length; e++) {
                                let t = l - a[e][1].length;
                                for (i = 0; i < t; i++) a[e][1].push([]);
                                for (a[e][1].push([o[e]]), i = 0; i < a[e][1].length; i++) {
                                    let t = g - a[e][1][i].length;
                                    for (var n = 0; n < t; n++) a[e][1][i].push(0)
                                }
                            }
                            a.sort((function (t, e) {
                                for (var i = 0; i < t[1].length; i++)
                                    for (var r = 0; r < t[1][i].length; r++) {
                                        if (t[1][i][r] > e[1][i][r]) return -1;
                                        if (t[1][i][r] < e[1][i][r]) return 1
                                    }
                                return 0
                            }));
                            let d = new Uint8Array(h);
                            for (e = 0; e < h; e++) d[e] = a[e][0], s.value.priority = e;
                            let u = this.graph.vertices[o[d[0]]].position,
                                c = this.graph.vertices[o[d[1]]].position,
                                p = this.graph.vertices[o[d[2]]].position,
                                f = u.relativeClockwise(c, s.position),
                                v = (u.relativeClockwise(p, s.position), -1 === f),
                                m = "@" === s.value.bracket.chirality ? -1 : 1,
                                b = r.parityOfPermutation(d) * m == 1 ? "R" : "S",
                                y = "down",
                                x = "up";
                            (v && "R" !== b || !v && "S" !== b) && (s.value.hydrogenDirection = "up", y = "up", x = "down"), s.value.hasHydrogen && (this.graph.getEdge(s.id, o[d[d.length - 1]]).wedge = y);
                            let S = new Array(o.length - 1),
                                A = s.value.rings.length > 1 && s.value.hasHydrogen,
                                C = s.value.hasHydrogen ? 1 : 0;
                            for (e = 0; e < d.length - C; e++) {
                                S[e] = new Uint32Array(2);
                                let t = this.graph.vertices[o[d[e]]];
                                S[e][0] += t.value.isStereoCenter ? 0 : 1e5, S[e][0] += this.areVerticesInSameRing(t, s) ? 0 : 1e4, S[e][0] += t.value.isHeteroAtom() ? 1e3 : 0, S[e][0] -= 0 === t.value.subtreeDepth ? 1e3 : 0, S[e][0] += 1e3 - t.value.subtreeDepth, S[e][1] = o[d[e]]
                            }
                            if (S.sort((function (t, e) {
                                return t[0] > e[0] ? -1 : t[0] < e[0] ? 1 : 0
                            })), !A) {
                                let t = S[0][1];
                                if (s.value.hasHydrogen) this.graph.getEdge(s.id, t).wedge = x;
                                else {
                                    let i = x;
                                    for (e = d.length - 1; e >= 0 && (i = i === y ? x : y, o[d[e]] !== t); e--);
                                    this.graph.getEdge(s.id, t).wedge = i
                                }
                            }
                            s.value.chirality = b
                        }
                    }
                    visitStereochemistry(t, e, i, r, n, s, o = 0) {
                        i[t] = 1;
                        let h = this.graph.vertices[t],
                            a = h.value.getAtomicNumber();
                        r.length <= s && r.push(Array());
                        for (var l = 0; l < this.graph.getEdge(t, e).weight; l++) r[s].push(1e3 * o + a);
                        let g = this.graph.vertices[t].neighbours;
                        for (l = 0; l < g.length; l++) 1 !== i[g[l]] && s < n - 1 && this.visitStereochemistry(g[l], t, i.slice(), r, n, s + 1, a);
                        if (s < n - 1) {
                            let e = 0;
                            for (l = 0; l < g.length; l++) e += this.graph.getEdge(t, g[l]).weight;
                            for (l = 0; l < h.value.getMaxBonds() - e; l++) r.length <= s + 1 && r.push(Array()), r[s + 1].push(1e3 * a + 1)
                        }
                    }
                    initPseudoElements() {
                        for (var t = 0; t < this.graph.vertices.length; t++) {
                            const i = this.graph.vertices[t],
                                r = i.neighbours;
                            let n = Array(r.length);
                            for (var e = 0; e < r.length; e++) n[e] = this.graph.vertices[r[e]];
                            if (i.getNeighbourCount() < 3 || i.value.rings.length > 0) continue;
                            if ("P" === i.value.element) continue;
                            if ("C" === i.value.element && 3 === n.length && "N" === n[0].value.element && "N" === n[1].value.element && "N" === n[2].value.element) continue;
                            let s = 0,
                                o = 0;
                            for (e = 0; e < n.length; e++) {
                                let t = n[e],
                                    i = t.value.element,
                                    r = t.getNeighbourCount();
                                "C" !== i && "H" !== i && 1 === r && s++, r > 1 && o++
                            }
                            if (o > 1 || s < 2) continue;
                            let h = null;
                            for (e = 0; e < n.length; e++) {
                                let t = n[e];
                                t.getNeighbourCount() > 1 && (h = t)
                            }
                            for (e = 0; e < n.length; e++) {
                                let t = n[e];
                                if (t.getNeighbourCount() > 1) continue;
                                t.value.isDrawn = !1;
                                let r = a.maxBonds[t.value.element] - t.value.bondCount,
                                    s = "";
                                t.value.bracket && (r = t.value.bracket.hcount, s = t.value.bracket.charge || 0), i.value.attachPseudoElement(t.value.element, h ? h.value.element : null, r, s)
                            }
                        }
                        for (t = 0; t < this.graph.vertices.length; t++) {
                            const i = this.graph.vertices[t],
                                r = i.value,
                                n = r.element;
                            if ("C" === n || "H" === n || !r.isDrawn) continue;
                            const s = i.neighbours;
                            let o = Array(s.length);
                            for (e = 0; e < s.length; e++) o[e] = this.graph.vertices[s[e]];
                            for (e = 0; e < o.length; e++) {
                                let t = o[e].value;
                                if (!t.hasAttachedPseudoElements || 2 !== t.getAttachedPseudoElementsCount()) continue;
                                const r = t.getAttachedPseudoElements();
                                r.hasOwnProperty("0O") && r.hasOwnProperty("3C") && (t.isDrawn = !1, i.value.attachPseudoElement("Ac", "", 0))
                            }
                        }
                    }
                }
            },
            826: t => {
                class e {
                    constructor(t, e, i = 1) {
                        this.id = null, this.sourceId = t, this.targetId = e, this.weight = i, this.bondType = "-", this.isPartOfAromaticRing = !1, this.center = !1, this.wedge = ""
                    }
                    setBondType(t) {
                        this.bondType = t, this.weight = e.bonds[t]
                    }
                    static get bonds() {
                        return {
                            "-": 1,
                            "/": 1,
                            "\\": 1,
                            "=": 2,
                            "#": 3,
                            $: 4
                        }
                    }
                }
                t.exports = e
            },
            707: (t, e, i) => {
                const r = i(474),
                    n = (i(614), i(843)),
                    s = i(826),
                    o = (i(421), i(427));
                class h {
                    constructor(t, e = !1) {
                        this.vertices = Array(), this.edges = Array(), this.vertexIdsToEdgeId = {}, this.isomeric = e, this._time = 0, this._init(t)
                    }
                    _init(t, e = 0, i = null, r = !1) {
                        let h = new o(t.atom.element ? t.atom.element : t.atom, t.bond);
                        h.branchBond = t.branchBond, h.ringbonds = t.ringbonds, h.bracket = t.atom.element ? t.atom : null, h.class = t.atom.class;
                        let a = new n(h),
                            l = this.vertices[i];
                        if (this.addVertex(a), null !== i) {
                            a.setParentVertexId(i), a.value.addNeighbouringElement(l.value.element), l.addChild(a.id), l.value.addNeighbouringElement(h.element), l.spanningTreeChildren.push(a.id);
                            let t = new s(i, a.id, 1),
                                e = null;
                            r ? (t.setBondType(a.value.branchBond || "-"), e = a.id, t.setBondType(a.value.branchBond || "-"), e = a.id) : (t.setBondType(l.value.bondType || "-"), e = l.id), this.addEdge(t)
                        }
                        let g = t.ringbondCount + 1;
                        h.bracket && (g += h.bracket.hcount);
                        let d = 0;
                        if (h.bracket && h.bracket.chirality) {
                            h.isStereoCenter = !0, d = h.bracket.hcount;
                            for (var u = 0; u < d; u++) this._init({
                                atom: "H",
                                isBracket: "false",
                                branches: Array(),
                                branchCount: 0,
                                ringbonds: Array(),
                                ringbondCount: !1,
                                next: null,
                                hasNext: !1,
                                bond: "-"
                            }, u, a.id, !0)
                        }
                        for (u = 0; u < t.branchCount; u++) this._init(t.branches[u], u + g, a.id, !0);
                        t.hasNext && this._init(t.next, t.branchCount + g, a.id)
                    }
                    clear() {
                        this.vertices = Array(), this.edges = Array(), this.vertexIdsToEdgeId = {}
                    }
                    addVertex(t) {
                        return t.id = this.vertices.length, this.vertices.push(t), t.id
                    }
                    addEdge(t) {
                        let e = this.vertices[t.sourceId],
                            i = this.vertices[t.targetId];
                        return t.id = this.edges.length, this.edges.push(t), this.vertexIdsToEdgeId[t.sourceId + "_" + t.targetId] = t.id, this.vertexIdsToEdgeId[t.targetId + "_" + t.sourceId] = t.id, t.isPartOfAromaticRing = e.value.isPartOfAromaticRing && i.value.isPartOfAromaticRing, e.value.bondCount += t.weight, i.value.bondCount += t.weight, e.edges.push(t.id), i.edges.push(t.id), t.id
                    }
                    getEdge(t, e) {
                        let i = this.vertexIdsToEdgeId[t + "_" + e];
                        return void 0 === i ? null : this.edges[i]
                    }
                    getEdges(t) {
                        let e = Array(),
                            i = this.vertices[t];
                        for (var r = 0; r < i.neighbours.length; r++) e.push(this.vertexIdsToEdgeId[t + "_" + i.neighbours[r]]);
                        return e
                    }
                    hasEdge(t, e) {
                        return void 0 !== this.vertexIdsToEdgeId[t + "_" + e]
                    }
                    getVertexList() {
                        let t = [this.vertices.length];
                        for (var e = 0; e < this.vertices.length; e++) t[e] = this.vertices[e].id;
                        return t
                    }
                    getEdgeList() {
                        let t = Array(this.edges.length);
                        for (var e = 0; e < this.edges.length; e++) t[e] = [this.edges[e].sourceId, this.edges[e].targetId];
                        return t
                    }
                    getAdjacencyMatrix() {
                        let t = this.vertices.length,
                            e = Array(t);
                        for (var i = 0; i < t; i++) e[i] = new Array(t), e[i].fill(0);
                        for (i = 0; i < this.edges.length; i++) {
                            let t = this.edges[i];
                            e[t.sourceId][t.targetId] = 1, e[t.targetId][t.sourceId] = 1
                        }
                        return e
                    }
                    getComponentsAdjacencyMatrix() {
                        let t = this.vertices.length,
                            e = Array(t),
                            i = this.getBridges();
                        for (var r = 0; r < t; r++) e[r] = new Array(t), e[r].fill(0);
                        for (r = 0; r < this.edges.length; r++) {
                            let t = this.edges[r];
                            e[t.sourceId][t.targetId] = 1, e[t.targetId][t.sourceId] = 1
                        }
                        for (r = 0; r < i.length; r++) e[i[r][0]][i[r][1]] = 0, e[i[r][1]][i[r][0]] = 0;
                        return e
                    }
                    getSubgraphAdjacencyMatrix(t) {
                        let e = t.length,
                            i = Array(e);
                        for (var r = 0; r < e; r++) {
                            i[r] = new Array(e), i[r].fill(0);
                            for (var n = 0; n < e; n++) r !== n && this.hasEdge(t[r], t[n]) && (i[r][n] = 1)
                        }
                        return i
                    }
                    getDistanceMatrix() {
                        let t = this.vertices.length,
                            e = this.getAdjacencyMatrix(),
                            i = Array(t);
                        for (var r = 0; r < t; r++) i[r] = Array(t), i[r].fill(1 / 0);
                        for (r = 0; r < t; r++)
                            for (var n = 0; n < t; n++) 1 === e[r][n] && (i[r][n] = 1);
                        for (var s = 0; s < t; s++)
                            for (r = 0; r < t; r++)
                                for (n = 0; n < t; n++) i[r][n] > i[r][s] + i[s][n] && (i[r][n] = i[r][s] + i[s][n]);
                        return i
                    }
                    getSubgraphDistanceMatrix(t) {
                        let e = t.length,
                            i = this.getSubgraphAdjacencyMatrix(t),
                            r = Array(e);
                        for (var n = 0; n < e; n++) r[n] = Array(e), r[n].fill(1 / 0);
                        for (n = 0; n < e; n++)
                            for (var s = 0; s < e; s++) 1 === i[n][s] && (r[n][s] = 1);
                        for (var o = 0; o < e; o++)
                            for (n = 0; n < e; n++)
                                for (s = 0; s < e; s++) r[n][s] > r[n][o] + r[o][s] && (r[n][s] = r[n][o] + r[o][s]);
                        return r
                    }
                    getAdjacencyList() {
                        let t = this.vertices.length,
                            e = Array(t);
                        for (var i = 0; i < t; i++) {
                            e[i] = [];
                            for (var r = 0; r < t; r++) i !== r && this.hasEdge(this.vertices[i].id, this.vertices[r].id) && e[i].push(r)
                        }
                        return e
                    }
                    getSubgraphAdjacencyList(t) {
                        let e = t.length,
                            i = Array(e);
                        for (var r = 0; r < e; r++) {
                            i[r] = Array();
                            for (var n = 0; n < e; n++) r !== n && this.hasEdge(t[r], t[n]) && i[r].push(n)
                        }
                        return i
                    }
                    getBridges() {
                        let t = this.vertices.length,
                            e = new Array(t),
                            i = new Array(t),
                            r = new Array(t),
                            n = new Array(t),
                            s = this.getAdjacencyList(),
                            o = Array();
                        e.fill(!1), n.fill(null), this._time = 0;
                        for (var h = 0; h < t; h++) e[h] || this._bridgeDfs(h, e, i, r, n, s, o);
                        return o
                    }
                    traverseBF(t, e) {
                        let i = this.vertices.length,
                            r = new Array(i);
                        r.fill(!1);
                        for (var n = [t]; n.length > 0;) {
                            let t = n.shift(),
                                i = this.vertices[t];
                            e(i);
                            for (var s = 0; s < i.neighbours.length; s++) {
                                let t = i.neighbours[s];
                                r[t] || (r[t] = !0, n.push(t))
                            }
                        }
                    }
                    getTreeDepth(t, e) {
                        if (null === t || null === e) return 0;
                        let i = this.vertices[t].getSpanningTreeNeighbours(e),
                            r = 0;
                        for (var n = 0; n < i.length; n++) {
                            let e = i[n],
                                s = this.getTreeDepth(e, t);
                            s > r && (r = s)
                        }
                        return r + 1
                    }
                    traverseTree(t, e, i, r = 999999, n = !1, s = 1, o = null) {
                        if (null === o && (o = new Uint8Array(this.vertices.length)), s > r + 1 || 1 === o[t]) return;
                        o[t] = 1;
                        let h = this.vertices[t],
                            a = h.getNeighbours(e);
                        (!n || s > 1) && i(h);
                        for (var l = 0; l < a.length; l++) this.traverseTree(a[l], t, i, r, n, s + 1, o)
                    }
                    kkLayout(t, e, i, n, s, o = .1, h = .1, a = 2e3, l = 50, g = 1e9) {
                        let d = s;
                        for (var u = t.length; u--;) var c = this.vertices[t[u]].neighbours.length;
                        let p = this.getSubgraphDistanceMatrix(t),
                            f = t.length,
                            v = r.polyCircumradius(500, f),
                            m = r.centralAngle(f),
                            b = 0,
                            y = new Float32Array(f),
                            x = new Float32Array(f),
                            S = Array(f);
                        for (u = f; u--;) {
                            let i = this.vertices[t[u]];
                            i.positioned ? (y[u] = i.position.x, x[u] = i.position.y) : (y[u] = e.x + Math.cos(b) * v, x[u] = e.y + Math.sin(b) * v), S[u] = i.positioned, b += m
                        }
                        let A = Array(f);
                        for (u = f; u--;)
                            for (A[u] = new Array(f), c = f; c--;) A[u][c] = s * p[u][c];
                        let C = Array(f);
                        for (u = f; u--;)
                            for (C[u] = Array(f), c = f; c--;) C[u][c] = d * Math.pow(p[u][c], -2);
                        let R, w, T, B, I, P, L, N = Array(f),
                            O = new Float32Array(f),
                            M = new Float32Array(f);
                        for (u = f; u--;) N[u] = Array(f);
                        for (u = f; u--;) {
                            R = y[u], w = x[u], T = 0, B = 0;
                            let t = f;
                            for (; t--;) u !== t && (I = y[t], P = x[t], L = 1 / Math.sqrt((R - I) * (R - I) + (w - P) * (w - P)), N[u][t] = [C[u][t] * (R - I - A[u][t] * (R - I) * L), C[u][t] * (w - P - A[u][t] * (w - P) * L)], N[t][u] = N[u][t], T += N[u][t][0], B += N[u][t][1]);
                            O[u] = T, M[u] = B
                        }
                        let k = function (t) {
                            return [O[t] * O[t] + M[t] * M[t], O[t], M[t]]
                        },
                            E = function () {
                                let t = 0,
                                    e = 0,
                                    i = 0,
                                    r = 0;
                                for (u = f; u--;) {
                                    let [n, s, o] = k(u);
                                    n > t && !1 === S[u] && (t = n, e = u, i = s, r = o)
                                }
                                return [e, t, i, r]
                            },
                            D = function (t, e, i) {
                                let r = 0,
                                    n = 0,
                                    s = 0,
                                    o = y[t],
                                    h = x[t],
                                    a = A[t],
                                    l = C[t];
                                for (u = f; u--;) {
                                    if (u === t) continue;
                                    let e = y[u],
                                        i = x[u],
                                        g = a[u],
                                        d = l[u],
                                        c = (o - e) * (o - e),
                                        p = 1 / Math.pow(c + (h - i) * (h - i), 1.5);
                                    r += d * (1 - g * (h - i) * (h - i) * p), n += d * (1 - g * c * p), s += d * (g * (o - e) * (h - i) * p)
                                }
                                0 === r && (r = .1), 0 === n && (n = .1), 0 === s && (s = .1);
                                let g = e / r + i / s;
                                g /= s / r - n / s;
                                let d = -(s * g + e) / r;
                                y[t] += d, x[t] += g;
                                let c, p, v, m, b, S = N[t];
                                for (e = 0, i = 0, o = y[t], h = x[t], u = f; u--;) t !== u && (c = y[u], p = x[u], v = S[u][0], m = S[u][1], b = 1 / Math.sqrt((o - c) * (o - c) + (h - p) * (h - p)), d = l[u] * (o - c - a[u] * (o - c) * b), g = l[u] * (h - p - a[u] * (h - p) * b), S[u] = [d, g], e += d, i += g, O[u] += d - v, M[u] += g - m);
                                O[t] = e, M[t] = i
                            },
                            F = 0,
                            z = 0,
                            H = 0,
                            V = 0,
                            W = 0,
                            U = 0;
                        for (; g > o && a > W;)
                            for (W++, [F, g, z, H] = E(), V = g, U = 0; V > h && l > U;) U++, D(F, z, H), [V, z, H] = k(F);
                        for (u = f; u--;) {
                            let e = t[u],
                                i = this.vertices[e];
                            i.position.x = y[u], i.position.y = x[u], i.positioned = !0, i.forcePositioned = !0
                        }
                    }
                    _bridgeDfs(t, e, i, r, n, s, o) {
                        e[t] = !0, i[t] = r[t] = ++this._time;
                        for (var h = 0; h < s[t].length; h++) {
                            let a = s[t][h];
                            e[a] ? a !== n[t] && (r[t] = Math.min(r[t], i[a])) : (n[a] = t, this._bridgeDfs(a, e, i, r, n, s, o), r[t] = Math.min(r[t], r[a]), r[a] > i[t] && o.push([t, a]))
                        }
                    }
                    static getConnectedComponents(t) {
                        let e = t.length,
                            i = new Array(e),
                            r = new Array;
                        i.fill(!1);
                        for (var n = 0; n < e; n++)
                            if (!i[n]) {
                                let e = Array();
                                i[n] = !0, e.push(n), h._ccGetDfs(n, i, t, e), e.length > 1 && r.push(e)
                            } return r
                    }
                    static getConnectedComponentCount(t) {
                        let e = t.length,
                            i = new Array(e),
                            r = 0;
                        i.fill(!1);
                        for (var n = 0; n < e; n++) i[n] || (i[n] = !0, r++, h._ccCountDfs(n, i, t));
                        return r
                    }
                    static _ccCountDfs(t, e, i) {
                        for (var r = 0; r < i[t].length; r++) i[t][r] && !e[r] && t !== r && (e[r] = !0, h._ccCountDfs(r, e, i))
                    }
                    static _ccGetDfs(t, e, i, r) {
                        for (var n = 0; n < i[t].length; n++) i[t][n] && !e[n] && t !== n && (e[n] = !0, r.push(n), h._ccGetDfs(n, e, i, r))
                    }
                }
                t.exports = h
            },
            929: (t, e, i) => {
                const r = i(614);
                class n {
                    constructor(t = new r(0, 0), e = new r(0, 0), i = null, n = null, s = !1, o = !1) {
                        this.from = t, this.to = e, this.elementFrom = i, this.elementTo = n, this.chiralFrom = s, this.chiralTo = o
                    }
                    clone() {
                        return new n(this.from.clone(), this.to.clone(), this.elementFrom, this.elementTo)
                    }
                    getLength() {
                        return Math.sqrt(Math.pow(this.to.x - this.from.x, 2) + Math.pow(this.to.y - this.from.y, 2))
                    }
                    getAngle() {
                        return r.subtract(this.getRightVector(), this.getLeftVector()).angle()
                    }
                    getRightVector() {
                        return this.from.x < this.to.x ? this.to : this.from
                    }
                    getLeftVector() {
                        return this.from.x < this.to.x ? this.from : this.to
                    }
                    getRightElement() {
                        return this.from.x < this.to.x ? this.elementTo : this.elementFrom
                    }
                    getLeftElement() {
                        return this.from.x < this.to.x ? this.elementFrom : this.elementTo
                    }
                    getRightChiral() {
                        return this.from.x < this.to.x ? this.chiralTo : this.chiralFrom
                    }
                    getLeftChiral() {
                        return this.from.x < this.to.x ? this.chiralFrom : this.chiralTo
                    }
                    setRightVector(t, e) {
                        return this.from.x < this.to.x ? (this.to.x = t, this.to.y = e) : (this.from.x = t, this.from.y = e), this
                    }
                    setLeftVector(t, e) {
                        return this.from.x < this.to.x ? (this.from.x = t, this.from.y = e) : (this.to.x = t, this.to.y = e), this
                    }
                    rotateToXAxis() {
                        let t = this.getLeftVector();
                        return this.setRightVector(t.x + this.getLength(), t.y), this
                    }
                    rotate(t) {
                        let e = this.getLeftVector(),
                            i = this.getRightVector(),
                            r = Math.sin(t),
                            n = Math.cos(t),
                            s = n * (i.x - e.x) - r * (i.y - e.y) + e.x,
                            o = r * (i.x - e.x) - n * (i.y - e.y) + e.y;
                        return this.setRightVector(s, o), this
                    }
                    shortenFrom(t) {
                        let e = r.subtract(this.to, this.from);
                        return e.normalize(), e.multiplyScalar(t), this.from.add(e), this
                    }
                    shortenTo(t) {
                        let e = r.subtract(this.from, this.to);
                        return e.normalize(), e.multiplyScalar(t), this.to.add(e), this
                    }
                    shortenRight(t) {
                        return this.from.x < this.to.x ? this.shortenTo(t) : this.shortenFrom(t), this
                    }
                    shortenLeft(t) {
                        return this.from.x < this.to.x ? this.shortenFrom(t) : this.shortenTo(t), this
                    }
                    shorten(t) {
                        let e = r.subtract(this.from, this.to);
                        return e.normalize(), e.multiplyScalar(t / 2), this.to.add(e), this.from.subtract(e), this
                    }
                }
                t.exports = n
            },
            474: t => {
                class e {
                    static round(t, e) {
                        return e = e || 1, Number(Math.round(t + "e" + e) + "e-" + e)
                    }
                    static meanAngle(t) {
                        let e = 0,
                            i = 0;
                        for (var r = 0; r < t.length; r++) e += Math.sin(t[r]), i += Math.cos(t[r]);
                        return Math.atan2(e / t.length, i / t.length)
                    }
                    static innerAngle(t) {
                        return e.toRad(180 * (t - 2) / t)
                    }
                    static polyCircumradius(t, e) {
                        return t / (2 * Math.sin(Math.PI / e))
                    }
                    static apothem(t, e) {
                        return t * Math.cos(Math.PI / e)
                    }
                    static apothemFromSideLength(t, i) {
                        let r = e.polyCircumradius(t, i);
                        return e.apothem(r, i)
                    }
                    static centralAngle(t) {
                        return e.toRad(360 / t)
                    }
                    static toDeg(t) {
                        return t * e.degFactor
                    }
                    static toRad(t) {
                        return t * e.radFactor
                    }
                    static parityOfPermutation(t) {
                        let e = new Uint8Array(t.length),
                            i = 0,
                            r = function (i, n = 0) {
                                return 1 === e[i] ? n : (n++, e[i] = 1, r(t[i], n))
                            };
                        for (var n = 0; n < t.length; n++) 1 !== e[n] && (i += 1 - r(n) % 2);
                        return i % 2 ? -1 : 1
                    }
                    static get radFactor() {
                        return Math.PI / 180
                    }
                    static get degFactor() {
                        return 180 / Math.PI
                    }
                    static get twoPI() {
                        return 2 * Math.PI
                    }
                }
                t.exports = e
            },
            207: t => {
                t.exports = class {
                    static extend() {
                        let t = this,
                            e = {},
                            i = !1,
                            r = 0,
                            n = arguments.length;
                        "[object Boolean]" === Object.prototype.toString.call(arguments[0]) && (i = arguments[0], r++);
                        let s = function (r) {
                            for (var n in r) Object.prototype.hasOwnProperty.call(r, n) && (i && "[object Object]" === Object.prototype.toString.call(r[n]) ? e[n] = t.extend(!0, e[n], r[n]) : e[n] = r[n])
                        };
                        for (; r < n; r++) s(arguments[r]);
                        return e
                    }
                }
            },
            19: t => {
                t.exports = function () {
                    "use strict";

                    function t(e, i, r, n) {
                        this.message = e, this.expected = i, this.found = r, this.location = n, this.name = "SyntaxError", "function" == typeof Error.captureStackTrace && Error.captureStackTrace(this, t)
                    }
                    return function (t, e) {
                        function i() {
                            this.constructor = t
                        }
                        i.prototype = e.prototype, t.prototype = new i
                    }(t, Error), t.buildMessage = function (t, e) {
                        var i = {
                            literal: function (t) {
                                return '"' + n(t.text) + '"'
                            },
                            class: function (t) {
                                var e, i = "";
                                for (e = 0; e < t.parts.length; e++) i += t.parts[e] instanceof Array ? s(t.parts[e][0]) + "-" + s(t.parts[e][1]) : s(t.parts[e]);
                                return "[" + (t.inverted ? "^" : "") + i + "]"
                            },
                            any: function (t) {
                                return "any character"
                            },
                            end: function (t) {
                                return "end of input"
                            },
                            other: function (t) {
                                return t.description
                            }
                        };

                        function r(t) {
                            return t.charCodeAt(0).toString(16).toUpperCase()
                        }

                        function n(t) {
                            return t.replace(/\\/g, "\\\\").replace(/"/g, '\\"').replace(/\0/g, "\\0").replace(/\t/g, "\\t").replace(/\n/g, "\\n").replace(/\r/g, "\\r").replace(/[\x00-\x0F]/g, (function (t) {
                                return "\\x0" + r(t)
                            })).replace(/[\x10-\x1F\x7F-\x9F]/g, (function (t) {
                                return "\\x" + r(t)
                            }))
                        }

                        function s(t) {
                            return t.replace(/\\/g, "\\\\").replace(/\]/g, "\\]").replace(/\^/g, "\\^").replace(/-/g, "\\-").replace(/\0/g, "\\0").replace(/\t/g, "\\t").replace(/\n/g, "\\n").replace(/\r/g, "\\r").replace(/[\x00-\x0F]/g, (function (t) {
                                return "\\x0" + r(t)
                            })).replace(/[\x10-\x1F\x7F-\x9F]/g, (function (t) {
                                return "\\x" + r(t)
                            }))
                        }
                        return "Expected " + function (t) {
                            var e, r, n, s = new Array(t.length);
                            for (e = 0; e < t.length; e++) s[e] = (n = t[e], i[n.type](n));
                            if (s.sort(), s.length > 0) {
                                for (e = 1, r = 1; e < s.length; e++) s[e - 1] !== s[e] && (s[r] = s[e], r++);
                                s.length = r
                            }
                            switch (s.length) {
                                case 1:
                                    return s[0];
                                case 2:
                                    return s[0] + " or " + s[1];
                                default:
                                    return s.slice(0, -1).join(", ") + ", or " + s[s.length - 1]
                            }
                        }(t) + " but " + function (t) {
                            return t ? '"' + n(t) + '"' : "end of input"
                        }(e) + " found."
                    }, {
                        SyntaxError: t,
                        parse: function (e, i) {
                            if (i = void 0 !== i ? i : {}, e.split("(").length - 1 != e.split(")").length - 1) throw ht("The number of opening parentheses does not match the number of closing parentheses.", 0);
                            var r, n, s, o, h = {},
                                a = {
                                    chain: at
                                },
                                l = at,
                                g = it("(", !1),
                                d = it(")", !1),
                                u = /^[\-=#$:\/\\.]/,
                                c = rt(["-", "=", "#", "$", ":", "/", "\\", "."], !1, !1),
                                p = it("[", !1),
                                f = it("se", !1),
                                v = it("as", !1),
                                m = it("]", !1),
                                b = it("B", !1),
                                y = it("r", !1),
                                x = it("C", !1),
                                S = it("l", !1),
                                A = /^[NOPSFI]/,
                                C = rt(["N", "O", "P", "S", "F", "I"], !1, !1),
                                R = /^[bcnops]/,
                                w = rt(["b", "c", "n", "o", "p", "s"], !1, !1),
                                T = it("*", !1),
                                B = /^[A-Z]/,
                                I = rt([
                                    ["A", "Z"]
                                ], !1, !1),
                                P = /^[a-z]/,
                                L = rt([
                                    ["a", "z"]
                                ], !1, !1),
                                N = it("%", !1),
                                O = /^[1-9]/,
                                M = rt([
                                    ["1", "9"]
                                ], !1, !1),
                                k = /^[0-9]/,
                                E = rt([
                                    ["0", "9"]
                                ], !1, !1),
                                D = it("@", !1),
                                F = it("TH", !1),
                                z = /^[12]/,
                                H = rt(["1", "2"], !1, !1),
                                V = it("AL", !1),
                                W = it("SP", !1),
                                U = /^[1-3]/,
                                q = rt([
                                    ["1", "3"]
                                ], !1, !1),
                                j = it("TB", !1),
                                _ = it("OH", !1),
                                G = it("+", !1),
                                X = it("-", !1),
                                K = it("H", !1),
                                Y = it(":", !1),
                                Z = /^[0]/,
                                $ = rt(["0"], !1, !1),
                                J = 0,
                                Q = [{
                                    line: 1,
                                    column: 1
                                }],
                                tt = 0,
                                et = [];
                            if ("startRule" in i) {
                                if (!(i.startRule in a)) throw new Error("Can't start parsing from rule \"" + i.startRule + '".');
                                l = a[i.startRule]
                            }

                            function it(t, e) {
                                return {
                                    type: "literal",
                                    text: t,
                                    ignoreCase: e
                                }
                            }

                            function rt(t, e, i) {
                                return {
                                    type: "class",
                                    parts: t,
                                    inverted: e,
                                    ignoreCase: i
                                }
                            }

                            function nt(t) {
                                var i, r = Q[t];
                                if (r) return r;
                                for (i = t - 1; !Q[i];) i--;
                                for (r = {
                                    line: (r = Q[i]).line,
                                    column: r.column
                                }; i < t;) 10 === e.charCodeAt(i) ? (r.line++, r.column = 1) : r.column++, i++;
                                return Q[t] = r, r
                            }

                            function st(t, e) {
                                var i = nt(t),
                                    r = nt(e);
                                return {
                                    start: {
                                        offset: t,
                                        line: i.line,
                                        column: i.column
                                    },
                                    end: {
                                        offset: e,
                                        line: r.line,
                                        column: r.column
                                    }
                                }
                            }

                            function ot(t) {
                                J < tt || (J > tt && (tt = J, et = []), et.push(t))
                            }

                            function ht(e, i) {
                                return new t(e, null, null, i)
                            }

                            function at() {
                                var t, i, r, n, s, o, a, l, g;
                                if (J, t = J, i = function () {
                                    var t;
                                    return J, t = function () {
                                        var t, i, r, n;
                                        return J, t = J, 66 === e.charCodeAt(J) ? (i = "B", J++) : (i = h, ot(b)), i !== h ? (114 === e.charCodeAt(J) ? (r = "r", J++) : (r = h, ot(y)), r === h && (r = null), r !== h ? t = i = [i, r] : (J = t, t = h)) : (J = t, t = h), t === h && (t = J, 67 === e.charCodeAt(J) ? (i = "C", J++) : (i = h, ot(x)), i !== h ? (108 === e.charCodeAt(J) ? (r = "l", J++) : (r = h, ot(S)), r === h && (r = null), r !== h ? t = i = [i, r] : (J = t, t = h)) : (J = t, t = h), t === h && (A.test(e.charAt(J)) ? (t = e.charAt(J), J++) : (t = h, ot(C)))), t !== h && (t = (n = t).length > 1 ? n.join("") : n), t
                                    }(), t === h && (t = dt()) === h && (t = function () {
                                        var t, i, r, n, s, o, a, l, g, d;
                                        return J, t = J, 91 === e.charCodeAt(J) ? (i = "[", J++) : (i = h, ot(p)), i !== h ? (r = function () {
                                            var t, i, r, n;
                                            return J, t = J, O.test(e.charAt(J)) ? (i = e.charAt(J), J++) : (i = h, ot(M)), i !== h ? (k.test(e.charAt(J)) ? (r = e.charAt(J), J++) : (r = h, ot(E)), r === h && (r = null), r !== h ? (k.test(e.charAt(J)) ? (n = e.charAt(J), J++) : (n = h, ot(E)), n === h && (n = null), n !== h ? t = i = [i, r, n] : (J = t, t = h)) : (J = t, t = h)) : (J = t, t = h), t !== h && (t = Number(t.join(""))), t
                                        }(), r === h && (r = null), r !== h ? ("se" === e.substr(J, 2) ? (n = "se", J += 2) : (n = h, ot(f)), n === h && ("as" === e.substr(J, 2) ? (n = "as", J += 2) : (n = h, ot(v)), n === h && (n = dt()) === h && (n = function () {
                                            var t, i, r;
                                            return J, t = J, B.test(e.charAt(J)) ? (i = e.charAt(J), J++) : (i = h, ot(I)), i !== h ? (P.test(e.charAt(J)) ? (r = e.charAt(J), J++) : (r = h, ot(L)), r === h && (r = null), r !== h ? t = i = [i, r] : (J = t, t = h)) : (J = t, t = h), t !== h && (t = t.join("")), t
                                        }(), n === h && (n = ut()))), n !== h ? (s = function () {
                                            var t, i, r, n, s, o, a;
                                            return J, t = J, 64 === e.charCodeAt(J) ? (i = "@", J++) : (i = h, ot(D)), i !== h ? (64 === e.charCodeAt(J) ? (r = "@", J++) : (r = h, ot(D)), r === h && (r = J, "TH" === e.substr(J, 2) ? (n = "TH", J += 2) : (n = h, ot(F)), n !== h ? (z.test(e.charAt(J)) ? (s = e.charAt(J), J++) : (s = h, ot(H)), s !== h ? r = n = [n, s] : (J = r, r = h)) : (J = r, r = h), r === h && (r = J, "AL" === e.substr(J, 2) ? (n = "AL", J += 2) : (n = h, ot(V)), n !== h ? (z.test(e.charAt(J)) ? (s = e.charAt(J), J++) : (s = h, ot(H)), s !== h ? r = n = [n, s] : (J = r, r = h)) : (J = r, r = h), r === h && (r = J, "SP" === e.substr(J, 2) ? (n = "SP", J += 2) : (n = h, ot(W)), n !== h ? (U.test(e.charAt(J)) ? (s = e.charAt(J), J++) : (s = h, ot(q)), s !== h ? r = n = [n, s] : (J = r, r = h)) : (J = r, r = h), r === h && (r = J, "TB" === e.substr(J, 2) ? (n = "TB", J += 2) : (n = h, ot(j)), n !== h ? (O.test(e.charAt(J)) ? (s = e.charAt(J), J++) : (s = h, ot(M)), s !== h ? (k.test(e.charAt(J)) ? (o = e.charAt(J), J++) : (o = h, ot(E)), o === h && (o = null), o !== h ? r = n = [n, s, o] : (J = r, r = h)) : (J = r, r = h)) : (J = r, r = h), r === h && (r = J, "OH" === e.substr(J, 2) ? (n = "OH", J += 2) : (n = h, ot(_)), n !== h ? (O.test(e.charAt(J)) ? (s = e.charAt(J), J++) : (s = h, ot(M)), s !== h ? (k.test(e.charAt(J)) ? (o = e.charAt(J), J++) : (o = h, ot(E)), o === h && (o = null), o !== h ? r = n = [n, s, o] : (J = r, r = h)) : (J = r, r = h)) : (J = r, r = h)))))), r === h && (r = null), r !== h ? t = i = [i, r] : (J = t, t = h)) : (J = t, t = h), t !== h && (t = (a = t)[1] ? "@" == a[1] ? "@@" : a[1].join("").replace(",", "") : "@"), t
                                        }(), s === h && (s = null), s !== h ? (o = function () {
                                            var t, i, r, n;
                                            return J, t = J, 72 === e.charCodeAt(J) ? (i = "H", J++) : (i = h, ot(K)), i !== h ? (k.test(e.charAt(J)) ? (r = e.charAt(J), J++) : (r = h, ot(E)), r === h && (r = null), r !== h ? t = i = [i, r] : (J = t, t = h)) : (J = t, t = h), t !== h && (t = (n = t)[1] ? Number(n[1]) : 1), t
                                        }(), o === h && (o = null), o !== h ? (a = function () {
                                            var t;
                                            return J, t = function () {
                                                var t, i, r, n, s, o;
                                                return J, t = J, 43 === e.charCodeAt(J) ? (i = "+", J++) : (i = h, ot(G)), i !== h ? (43 === e.charCodeAt(J) ? (r = "+", J++) : (r = h, ot(G)), r === h && (r = J, O.test(e.charAt(J)) ? (n = e.charAt(J), J++) : (n = h, ot(M)), n !== h ? (k.test(e.charAt(J)) ? (s = e.charAt(J), J++) : (s = h, ot(E)), s === h && (s = null), s !== h ? r = n = [n, s] : (J = r, r = h)) : (J = r, r = h)), r === h && (r = null), r !== h ? t = i = [i, r] : (J = t, t = h)) : (J = t, t = h), t !== h && (t = (o = t)[1] ? "+" != o[1] ? Number(o[1].join("")) : 2 : 1), t
                                            }(), t === h && (t = function () {
                                                var t, i, r, n, s, o;
                                                return J, t = J, 45 === e.charCodeAt(J) ? (i = "-", J++) : (i = h, ot(X)), i !== h ? (45 === e.charCodeAt(J) ? (r = "-", J++) : (r = h, ot(X)), r === h && (r = J, O.test(e.charAt(J)) ? (n = e.charAt(J), J++) : (n = h, ot(M)), n !== h ? (k.test(e.charAt(J)) ? (s = e.charAt(J), J++) : (s = h, ot(E)), s === h && (s = null), s !== h ? r = n = [n, s] : (J = r, r = h)) : (J = r, r = h)), r === h && (r = null), r !== h ? t = i = [i, r] : (J = t, t = h)) : (J = t, t = h), t !== h && (t = (o = t)[1] ? "-" != o[1] ? -Number(o[1].join("")) : -2 : -1), t
                                            }()), t
                                        }(), a === h && (a = null), a !== h ? (l = function () {
                                            var t, i, r, n, s, o, a;
                                            if (J, t = J, 58 === e.charCodeAt(J) ? (i = ":", J++) : (i = h, ot(Y)), i !== h) {
                                                if (r = J, O.test(e.charAt(J)) ? (n = e.charAt(J), J++) : (n = h, ot(M)), n !== h) {
                                                    for (s = [], k.test(e.charAt(J)) ? (o = e.charAt(J), J++) : (o = h, ot(E)); o !== h;) s.push(o), k.test(e.charAt(J)) ? (o = e.charAt(J), J++) : (o = h, ot(E));
                                                    s !== h ? r = n = [n, s] : (J = r, r = h)
                                                } else J = r, r = h;
                                                r === h && (Z.test(e.charAt(J)) ? (r = e.charAt(J), J++) : (r = h, ot($))), r !== h ? t = i = [i, r] : (J = t, t = h)
                                            } else J = t, t = h;
                                            return t !== h && (a = t, t = Number(a[1][0] + a[1][1].join(""))), t
                                        }(), l === h && (l = null), l !== h ? (93 === e.charCodeAt(J) ? (g = "]", J++) : (g = h, ot(m)), g !== h ? t = i = [i, r, n, s, o, a, l, g] : (J = t, t = h)) : (J = t, t = h)) : (J = t, t = h)) : (J = t, t = h)) : (J = t, t = h)) : (J = t, t = h)) : (J = t, t = h)) : (J = t, t = h), t !== h && (t = {
                                            isotope: (d = t)[1],
                                            element: d[2],
                                            chirality: d[3],
                                            hcount: d[4],
                                            charge: d[5],
                                            class: d[6]
                                        }), t
                                    }(), t === h && (t = ut())), t
                                }(), i !== h) {
                                    for (r = [], n = lt(); n !== h;) r.push(n), n = lt();
                                    if (r !== h) {
                                        for (n = [], s = J, (o = gt()) === h && (o = null), o !== h && (a = ct()) !== h ? s = o = [o, a] : (J = s, s = h); s !== h;) n.push(s), s = J, (o = gt()) === h && (o = null), o !== h && (a = ct()) !== h ? s = o = [o, a] : (J = s, s = h);
                                        if (n !== h) {
                                            for (s = [], o = lt(); o !== h;) s.push(o), o = lt();
                                            if (s !== h)
                                                if ((o = gt()) === h && (o = null), o !== h)
                                                    if ((a = at()) === h && (a = null), a !== h) {
                                                        for (l = [], g = lt(); g !== h;) l.push(g), g = lt();
                                                        l !== h ? t = i = [i, r, n, s, o, a, l] : (J = t, t = h)
                                                    } else J = t, t = h;
                                                else J = t, t = h;
                                            else J = t, t = h
                                        } else J = t, t = h
                                    } else J = t, t = h
                                } else J = t, t = h;
                                return t !== h && (t = function (t) {
                                    for (var e = [], i = [], r = 0; r < t[1].length; r++) e.push(t[1][r]);
                                    for (r = 0; r < t[2].length; r++) {
                                        var n = t[2][r][0] ? t[2][r][0] : "-";
                                        i.push({
                                            bond: n,
                                            id: t[2][r][1]
                                        })
                                    }
                                    for (r = 0; r < t[3].length; r++) e.push(t[3][r]);
                                    for (r = 0; r < t[6].length; r++) e.push(t[6][r]);
                                    return {
                                        atom: t[0],
                                        isBracket: !!t[0].element,
                                        branches: e,
                                        branchCount: e.length,
                                        ringbonds: i,
                                        ringbondCount: i.length,
                                        bond: t[4] ? t[4] : "-",
                                        next: t[5],
                                        hasNext: !!t[5]
                                    }
                                }(t)), t
                            }

                            function lt() {
                                var t, i, r, n, s, o, a;
                                return J, t = J, 40 === e.charCodeAt(J) ? (i = "(", J++) : (i = h, ot(g)), i !== h ? ((r = gt()) === h && (r = null), r !== h && (n = at()) !== h ? (41 === e.charCodeAt(J) ? (s = ")", J++) : (s = h, ot(d)), s !== h ? t = i = [i, r, n, s] : (J = t, t = h)) : (J = t, t = h)) : (J = t, t = h), t !== h && (a = (o = t)[1] ? o[1] : "-", o[2].branchBond = a, t = o[2]), t
                            }

                            function gt() {
                                var t;
                                if (J, u.test(e.charAt(J))) {
                                    if ((t = e.charAt(J)) === e.charAt(J + 1)) throw t = h, ht("The parser encountered a bond repetition.", J + 1);
                                    J++
                                } else t = h, ot(c);
                                return t
                            }

                            function dt() {
                                var t;
                                return J, R.test(e.charAt(J)) ? (t = e.charAt(J), J++) : (t = h, ot(w)), t
                            }

                            function ut() {
                                var t;
                                return J, 42 === e.charCodeAt(J) ? (t = "*", J++) : (t = h, ot(T)), t
                            }

                            function ct() {
                                var t, i, r, n, s;
                                return J, t = J, 37 === e.charCodeAt(J) ? (i = "%", J++) : (i = h, ot(N)), i !== h ? (O.test(e.charAt(J)) ? (r = e.charAt(J), J++) : (r = h, ot(M)), r !== h ? (k.test(e.charAt(J)) ? (n = e.charAt(J), J++) : (n = h, ot(E)), n !== h ? t = i = [i, r, n] : (J = t, t = h)) : (J = t, t = h)) : (J = t, t = h), t === h && (k.test(e.charAt(J)) ? (t = e.charAt(J), J++) : (t = h, ot(E))), t !== h && (t = 1 == (s = t).length ? Number(s) : Number(s.join("").replace("%", ""))), t
                            }
                            if ((r = l()) !== h && J === e.length) return r;
                            throw r !== h && J < e.length && ot({
                                type: "end"
                            }), n = et, s = tt < e.length ? e.charAt(tt) : null, o = tt < e.length ? st(tt, tt + 1) : st(tt, tt), new t(t.buildMessage(n, s), n, s, o)
                        }
                    }
                }()
            },
            421: (t, e, i) => {
                const r = i(348),
                    n = i(614),
                    s = (i(843), i(333));
                class o {
                    constructor(t) {
                        this.id = null, this.members = t, this.edges = [], this.insiders = [], this.neighbours = [], this.positioned = !1, this.center = new n(0, 0), this.rings = [], this.isBridged = !1, this.isPartOfBridged = !1, this.isSpiro = !1, this.isFused = !1, this.centralAngle = 0, this.canFlip = !0
                    }
                    clone() {
                        let t = new o(this.members);
                        return t.id = this.id, t.insiders = r.clone(this.insiders), t.neighbours = r.clone(this.neighbours), t.positioned = this.positioned, t.center = this.center.clone(), t.rings = r.clone(this.rings), t.isBridged = this.isBridged, t.isPartOfBridged = this.isPartOfBridged, t.isSpiro = this.isSpiro, t.isFused = this.isFused, t.centralAngle = this.centralAngle, t.canFlip = this.canFlip, t
                    }
                    getSize() {
                        return this.members.length
                    }
                    getPolygon(t) {
                        let e = [];
                        for (let i = 0; i < this.members.length; i++) e.push(t[this.members[i]].position);
                        return e
                    }
                    getAngle() {
                        return Math.PI - this.centralAngle
                    }
                    eachMember(t, e, i, r) {
                        let n = i = i || 0 === i ? i : this.members[0],
                            s = 0;
                        for (; null != n && s < 100;) {
                            let o = n;
                            e(o), n = t[n].getNextInRing(t, this.id, r), r = o, n == i && (n = null), s++
                        }
                    }
                    getOrderedNeighbours(t) {
                        let e = Array(this.neighbours.length);
                        for (let i = 0; i < this.neighbours.length; i++) {
                            let r = s.getVertices(t, this.id, this.neighbours[i]);
                            e[i] = {
                                n: r.length,
                                neighbour: this.neighbours[i]
                            }
                        }
                        return e.sort((function (t, e) {
                            return e.n - t.n
                        })), e
                    }
                    isBenzeneLike(t) {
                        let e = this.getDoubleBondCount(t),
                            i = this.members.length;
                        return 3 === e && 6 === i || 2 === e && 5 === i
                    }
                    getDoubleBondCount(t) {
                        let e = 0;
                        for (let i = 0; i < this.members.length; i++) {
                            let r = t[this.members[i]].value;
                            "=" !== r.bondType && "=" !== r.branchBond || e++
                        }
                        return e
                    }
                    contains(t) {
                        for (let e = 0; e < this.members.length; e++)
                            if (this.members[e] == t) return !0;
                        return !1
                    }
                }
                t.exports = o
            },
            333: (t, e, i) => {
                i(843), i(421), t.exports = class {
                    constructor(t, e) {
                        this.id = null, this.firstRingId = t.id, this.secondRingId = e.id, this.vertices = new Set;
                        for (var i = 0; i < t.members.length; i++) {
                            let r = t.members[i];
                            for (let t = 0; t < e.members.length; t++) r === e.members[t] && this.addVertex(r)
                        }
                    }
                    addVertex(t) {
                        this.vertices.add(t)
                    }
                    updateOther(t, e) {
                        this.firstRingId === e ? this.secondRingId = t : this.firstRingId = t
                    }
                    containsRing(t) {
                        return this.firstRingId === t || this.secondRingId === t
                    }
                    isBridge(t) {
                        if (this.vertices.size > 2) return !0;
                        for (let e of this.vertices)
                            if (t[e].value.rings.length > 2) return !0;
                        return !1
                    }
                    static isBridge(t, e, i, r) {
                        let n = null;
                        for (let s = 0; s < t.length; s++)
                            if (n = t[s], n.firstRingId === i && n.secondRingId === r || n.firstRingId === r && n.secondRingId === i) return n.isBridge(e);
                        return !1
                    }
                    static getNeighbours(t, e) {
                        let i = [];
                        for (let r = 0; r < t.length; r++) {
                            let n = t[r];
                            n.firstRingId === e ? i.push(n.secondRingId) : n.secondRingId === e && i.push(n.firstRingId)
                        }
                        return i
                    }
                    static getVertices(t, e, i) {
                        for (let r = 0; r < t.length; r++) {
                            let n = t[r];
                            if (n.firstRingId === e && n.secondRingId === i || n.firstRingId === i && n.secondRingId === e) return [...n.vertices]
                        }
                    }
                }
            },
            473: (t, e, i) => {
                const r = i(707);
                class n {
                    static getRings(t, e = !1) {
                        let i = t.getComponentsAdjacencyMatrix();
                        if (0 === i.length) return null;
                        let s = r.getConnectedComponents(i),
                            o = Array();
                        for (var h = 0; h < s.length; h++) {
                            let i = s[h],
                                r = t.getSubgraphAdjacencyMatrix([...i]),
                                g = new Uint16Array(r.length),
                                d = new Uint16Array(r.length);
                            for (var a = 0; a < r.length; a++) {
                                d[a] = 0, g[a] = 0;
                                for (var l = 0; l < r[a].length; l++) g[a] += r[a][l]
                            }
                            let u = 0;
                            for (a = 0; a < r.length; a++)
                                for (l = a + 1; l < r.length; l++) u += r[a][l];
                            let c = u - r.length + 1,
                                p = !0;
                            for (a = 0; a < g.length; a++) 3 !== g[a] && (p = !1);
                            if (p && (c = 2 + u - r.length), 1 === c) {
                                o.push([...i]);
                                continue
                            }
                            e && (c = 999);
                            let {
                                d: f,
                                pe: v,
                                pe_prime: m
                            } = n.getPathIncludedDistanceMatrices(r), b = n.getRingCandidates(f, v, m), y = n.getSSSR(b, f, r, v, m, g, d, c);
                            for (a = 0; a < y.length; a++) {
                                let t = Array(y[a].size),
                                    e = 0;
                                for (let r of y[a]) t[e++] = i[r];
                                o.push(t)
                            }
                        }
                        return o
                    }
                    static matrixToString(t) {
                        let e = "";
                        for (var i = 0; i < t.length; i++) {
                            for (var r = 0; r < t[i].length; r++) e += t[i][r] + " ";
                            e += "\n"
                        }
                        return e
                    }
                    static getPathIncludedDistanceMatrices(t) {
                        let e = t.length,
                            i = Array(e),
                            r = Array(e),
                            n = Array(e);
                        for (var s = 0, o = 0, h = 0, a = e; a--;) {
                            i[a] = Array(e), r[a] = Array(e), n[a] = Array(e);
                            for (var l = e; l--;) i[a][l] = a === l || 1 === t[a][l] ? t[a][l] : Number.POSITIVE_INFINITY, 1 === i[a][l] ? r[a][l] = [
                                [
                                    [a, l]
                                ]
                            ] : r[a][l] = Array(), n[a][l] = Array()
                        }
                        for (var g = e; g--;)
                            for (a = e; a--;)
                                for (l = e; l--;) {
                                    const t = i[a][l],
                                        e = i[a][g] + i[g][l];
                                    if (t > e) {
                                        if (t === e + 1)
                                            for (n[a][l] = [r[a][l].length], s = r[a][l].length; s--;)
                                                for (n[a][l][s] = [r[a][l][s].length], o = r[a][l][s].length; o--;)
                                                    for (n[a][l][s][o] = [r[a][l][s][o].length], h = r[a][l][s][o].length; h--;) n[a][l][s][o][h] = [r[a][l][s][o][0], r[a][l][s][o][1]];
                                        else n[a][l] = Array();
                                        for (i[a][l] = e, r[a][l] = [
                                            []
                                        ], s = r[a][g][0].length; s--;) r[a][l][0].push(r[a][g][0][s]);
                                        for (s = r[g][l][0].length; s--;) r[a][l][0].push(r[g][l][0][s])
                                    } else if (t === e) {
                                        if (r[a][g].length && r[g][l].length)
                                            if (r[a][l].length) {
                                                let t = Array();
                                                for (s = r[a][g][0].length; s--;) t.push(r[a][g][0][s]);
                                                for (s = r[g][l][0].length; s--;) t.push(r[g][l][0][s]);
                                                r[a][l].push(t)
                                            } else {
                                                let t = Array();
                                                for (s = r[a][g][0].length; s--;) t.push(r[a][g][0][s]);
                                                for (s = r[g][l][0].length; s--;) t.push(r[g][l][0][s]);
                                                r[a][l][0] = t
                                            }
                                    } else if (t === e - 1)
                                        if (n[a][l].length) {
                                            let t = Array();
                                            for (s = r[a][g][0].length; s--;) t.push(r[a][g][0][s]);
                                            for (s = r[g][l][0].length; s--;) t.push(r[g][l][0][s]);
                                            n[a][l].push(t)
                                        } else {
                                            let t = Array();
                                            for (s = r[a][g][0].length; s--;) t.push(r[a][g][0][s]);
                                            for (s = r[g][l][0].length; s--;) t.push(r[g][l][0][s]);
                                            n[a][l][0] = t
                                        }
                                }
                        return {
                            d: i,
                            pe: r,
                            pe_prime: n
                        }
                    }
                    static getRingCandidates(t, e, i) {
                        let r = t.length,
                            n = Array(),
                            s = 0;
                        for (let o = 0; o < r; o++)
                            for (let h = 0; h < r; h++) 0 === t[o][h] || 1 === e[o][h].length && 0 === i[o][h] || (s = 0 !== i[o][h].length ? 2 * (t[o][h] + .5) : 2 * t[o][h], s !== 1 / 0 && n.push([s, e[o][h], i[o][h]]));
                        return n.sort((function (t, e) {
                            return t[0] - e[0]
                        })), n
                    }
                    static getSSSR(t, e, i, r, s, o, h, a) {
                        let l = Array(),
                            g = Array();
                        for (let e = 0; e < t.length; e++)
                            if (t[e][0] % 2 != 0)
                                for (let r = 0; r < t[e][2].length; r++) {
                                    let s = t[e][1][0].concat(t[e][2][r]);
                                    for (var d = 0; d < s.length; d++) s[d][0].constructor === Array && (s[d] = s[d][0]);
                                    let u = n.bondsToAtoms(s);
                                    if (n.getBondCount(u, i) !== u.size || n.pathSetsContain(l, u, s, g, o, h) || (l.push(u), g = g.concat(s)), l.length > a) return l
                                } else
                                for (let r = 0; r < t[e][1].length - 1; r++) {
                                    let s = t[e][1][r].concat(t[e][1][r + 1]);
                                    for (d = 0; d < s.length; d++) s[d][0].constructor === Array && (s[d] = s[d][0]);
                                    let u = n.bondsToAtoms(s);
                                    if (n.getBondCount(u, i) !== u.size || n.pathSetsContain(l, u, s, g, o, h) || (l.push(u), g = g.concat(s)), l.length > a) return l
                                }
                        return l
                    }
                    static getEdgeCount(t) {
                        let e = 0,
                            i = t.length;
                        for (var r = i - 1; r--;)
                            for (var n = i; n--;) 1 === t[r][n] && e++;
                        return e
                    }
                    static getEdgeList(t) {
                        let e = t.length,
                            i = Array();
                        for (var r = e - 1; r--;)
                            for (var n = e; n--;) 1 === t[r][n] && i.push([r, n]);
                        return i
                    }
                    static bondsToAtoms(t) {
                        let e = new Set;
                        for (var i = t.length; i--;) e.add(t[i][0]), e.add(t[i][1]);
                        return e
                    }
                    static getBondCount(t, e) {
                        let i = 0;
                        for (let r of t)
                            for (let n of t) r !== n && (i += e[r][n]);
                        return i / 2
                    }
                    static pathSetsContain(t, e, i, r, s, o) {
                        for (var h = t.length; h--;) {
                            if (n.isSupersetOf(e, t[h])) return !0;
                            if (t[h].size === e.size && n.areSetsEqual(t[h], e)) return !0
                        }
                        let a = 0,
                            l = !1;
                        for (h = i.length; h--;)
                            for (var g = r.length; g--;)(i[h][0] === r[g][0] && i[h][1] === r[g][1] || i[h][1] === r[g][0] && i[h][0] === r[g][1]) && a++, a === i.length && (l = !0);
                        let d = !1;
                        if (l)
                            for (let t of e)
                                if (o[t] < s[t]) {
                                    d = !0;
                                    break
                                } if (l && !d) return !0;
                        for (let t of e) o[t]++;
                        return !1
                    }
                    static areSetsEqual(t, e) {
                        if (t.size !== e.size) return !1;
                        for (let i of t)
                            if (!e.has(i)) return !1;
                        return !0
                    }
                    static isSupersetOf(t, e) {
                        for (var i of e)
                            if (!t.has(i)) return !1;
                        return !0
                    }
                }
                t.exports = n
            },
            654: t => {
                t.exports = class {
                    constructor(t, e) {
                        this.colors = t, this.theme = this.colors[e]
                    }
                    getColor(t) {
                        return t && (t = t.toUpperCase()) in this.theme ? this.theme[t] : this.theme.C
                    }
                    setTheme(t) {
                        this.colors.hasOwnProperty(t) && (this.theme = this.colors[t])
                    }
                }
            },
            537: t => {
                t.exports = {
                    getChargeText: function (t) {
                        return 1 === t ? "+" : 2 === t ? "2+" : -1 === t ? "-" : -2 === t ? "2-" : ""
                    }
                }
            },
            614: t => {
                class e {
                    constructor(t, e) {
                        0 == arguments.length ? (this.x = 0, this.y = 0) : 1 == arguments.length ? (this.x = t.x, this.y = t.y) : (this.x = t, this.y = e)
                    }
                    clone() {
                        return new e(this.x, this.y)
                    }
                    toString() {
                        return "(" + this.x + "," + this.y + ")"
                    }
                    add(t) {
                        return this.x += t.x, this.y += t.y, this
                    }
                    subtract(t) {
                        return this.x -= t.x, this.y -= t.y, this
                    }
                    divide(t) {
                        return this.x /= t, this.y /= t, this
                    }
                    multiply(t) {
                        return this.x *= t.x, this.y *= t.y, this
                    }
                    multiplyScalar(t) {
                        return this.x *= t, this.y *= t, this
                    }
                    invert() {
                        return this.x = -this.x, this.y = -this.y, this
                    }
                    angle() {
                        return Math.atan2(this.y, this.x)
                    }
                    distance(t) {
                        return Math.sqrt((t.x - this.x) * (t.x - this.x) + (t.y - this.y) * (t.y - this.y))
                    }
                    distanceSq(t) {
                        return (t.x - this.x) * (t.x - this.x) + (t.y - this.y) * (t.y - this.y)
                    }
                    clockwise(t) {
                        let e = this.y * t.x,
                            i = this.x * t.y;
                        return e > i ? -1 : e === i ? 0 : 1
                    }
                    relativeClockwise(t, e) {
                        let i = (this.y - t.y) * (e.x - t.x),
                            r = (this.x - t.x) * (e.y - t.y);
                        return i > r ? -1 : i === r ? 0 : 1
                    }
                    rotate(t) {
                        let i = new e(0, 0),
                            r = Math.cos(t),
                            n = Math.sin(t);
                        return i.x = this.x * r - this.y * n, i.y = this.x * n + this.y * r, this.x = i.x, this.y = i.y, this
                    }
                    rotateAround(t, e) {
                        let i = Math.sin(t),
                            r = Math.cos(t);
                        this.x -= e.x, this.y -= e.y;
                        let n = this.x * r - this.y * i,
                            s = this.x * i + this.y * r;
                        return this.x = n + e.x, this.y = s + e.y, this
                    }
                    rotateTo(t, i, r = 0) {
                        this.x += .001, this.y -= .001;
                        let n = e.subtract(this, i),
                            s = e.subtract(t, i),
                            o = e.angle(s, n);
                        return this.rotateAround(o + r, i), this
                    }
                    rotateAwayFrom(t, e, i) {
                        this.rotateAround(i, e);
                        let r = this.distanceSq(t);
                        this.rotateAround(-2 * i, e), this.distanceSq(t) < r && this.rotateAround(2 * i, e)
                    }
                    getRotateAwayFromAngle(t, e, i) {
                        let r = this.clone();
                        r.rotateAround(i, e);
                        let n = r.distanceSq(t);
                        return r.rotateAround(-2 * i, e), r.distanceSq(t) < n ? i : -i
                    }
                    getRotateTowardsAngle(t, e, i) {
                        let r = this.clone();
                        r.rotateAround(i, e);
                        let n = r.distanceSq(t);
                        return r.rotateAround(-2 * i, e), r.distanceSq(t) > n ? i : -i
                    }
                    getRotateToAngle(t, i) {
                        let r = e.subtract(this, i),
                            n = e.subtract(t, i),
                            s = e.angle(n, r);
                        return Number.isNaN(s) ? 0 : s
                    }
                    isInPolygon(t) {
                        let e = !1;
                        for (let i = 0, r = t.length - 1; i < t.length; r = i++) t[i].y > this.y != t[r].y > this.y && this.x < (t[r].x - t[i].x) * (this.y - t[i].y) / (t[r].y - t[i].y) + t[i].x && (e = !e);
                        return e
                    }
                    length() {
                        return Math.sqrt(this.x * this.x + this.y * this.y)
                    }
                    lengthSq() {
                        return this.x * this.x + this.y * this.y
                    }
                    normalize() {
                        return this.divide(this.length()), this
                    }
                    normalized() {
                        return e.divideScalar(this, this.length())
                    }
                    whichSide(t, e) {
                        return (this.x - t.x) * (e.y - t.y) - (this.y - t.y) * (e.x - t.x)
                    }
                    sameSideAs(t, e, i) {
                        let r = this.whichSide(t, e),
                            n = i.whichSide(t, e);
                        return r < 0 && n < 0 || 0 == r && 0 == n || r > 0 && n > 0
                    }
                    static add(t, i) {
                        return new e(t.x + i.x, t.y + i.y)
                    }
                    static subtract(t, i) {
                        return new e(t.x - i.x, t.y - i.y)
                    }
                    static multiply(t, i) {
                        return new e(t.x * i.x, t.y * i.y)
                    }
                    static multiplyScalar(t, i) {
                        return new e(t.x, t.y).multiplyScalar(i)
                    }
                    static midpoint(t, i) {
                        return new e((t.x + i.x) / 2, (t.y + i.y) / 2)
                    }
                    static normals(t, i) {
                        let r = e.subtract(i, t);
                        return [new e(-r.y, r.x), new e(r.y, -r.x)]
                    }
                    static units(t, i) {
                        let r = e.subtract(i, t);
                        return [new e(-r.y, r.x).normalize(), new e(r.y, -r.x).normalize()]
                    }
                    static divide(t, i) {
                        return new e(t.x / i.x, t.y / i.y)
                    }
                    static divideScalar(t, i) {
                        return new e(t.x / i, t.y / i)
                    }
                    static dot(t, e) {
                        return t.x * e.x + t.y * e.y
                    }
                    static angle(t, i) {
                        let r = e.dot(t, i);
                        return Math.acos(r / (t.length() * i.length()))
                    }
                    static threePointangle(t, i, r) {
                        let n = e.subtract(i, t),
                            s = e.subtract(r, i),
                            o = t.distance(i),
                            h = i.distance(r);
                        return Math.acos(e.dot(n, s) / (o * h))
                    }
                    static scalarProjection(t, i) {
                        let r = i.normalized();
                        return e.dot(t, r)
                    }
                    static averageDirection(t) {
                        let i = new e(0, 0);
                        for (var r = 0; r < t.length; r++) {
                            let e = t[r];
                            i.add(e)
                        }
                        return i.normalize()
                    }
                }
                t.exports = e
            },
            843: (t, e, i) => {
                const r = i(474),
                    n = i(348),
                    s = i(614);
                i(427);
                class o {
                    constructor(t, e = 0, i = 0) {
                        this.id = null, this.value = t, this.position = new s(e || 0, i || 0), this.previousPosition = new s(0, 0), this.parentVertexId = null, this.children = Array(), this.spanningTreeChildren = Array(), this.edges = Array(), this.positioned = !1, this.angle = null, this.dir = 1, this.neighbourCount = 0, this.neighbours = Array(), this.neighbouringElements = Array(), this.forcePositioned = !1
                    }
                    setPosition(t, e) {
                        this.position.x = t, this.position.y = e
                    }
                    setPositionFromVector(t) {
                        this.position.x = t.x, this.position.y = t.y
                    }
                    addChild(t) {
                        this.children.push(t), this.neighbours.push(t), this.neighbourCount++
                    }
                    addRingbondChild(t, e) {
                        if (this.children.push(t), this.value.bracket) {
                            let i = 1;
                            0 === this.id && 0 === this.value.bracket.hcount && (i = 0), 1 === this.value.bracket.hcount && 0 === e && (i = 2), 1 === this.value.bracket.hcount && 1 === e && (i = this.neighbours.length < 3 ? 2 : 3), null === this.value.bracket.hcount && 0 === e && (i = 1), null === this.value.bracket.hcount && 1 === e && (i = this.neighbours.length < 3 ? 1 : 2), this.neighbours.splice(i, 0, t)
                        } else this.neighbours.push(t);
                        this.neighbourCount++
                    }
                    setParentVertexId(t) {
                        this.neighbourCount++, this.parentVertexId = t, this.neighbours.push(t)
                    }
                    isTerminal() {
                        return !!this.value.hasAttachedPseudoElements || null === this.parentVertexId && this.children.length < 2 || 0 === this.children.length
                    }
                    clone() {
                        let t = new o(this.value, this.position.x, this.position.y);
                        return t.id = this.id, t.previousPosition = new s(this.previousPosition.x, this.previousPosition.y), t.parentVertexId = this.parentVertexId, t.children = n.clone(this.children), t.spanningTreeChildren = n.clone(this.spanningTreeChildren), t.edges = n.clone(this.edges), t.positioned = this.positioned, t.angle = this.angle, t.forcePositioned = this.forcePositioned, t
                    }
                    equals(t) {
                        return this.id === t.id
                    }
                    getAngle(t = null, e = !1) {
                        let i = null;
                        return i = t ? s.subtract(this.position, t) : s.subtract(this.position, this.previousPosition), e ? r.toDeg(i.angle()) : i.angle()
                    }
                    getTextDirection(t) {
                        let e = this.getDrawnNeighbours(t),
                            i = Array();
                        if (1 === t.length) return "right";
                        for (let r = 0; r < e.length; r++) i.push(this.getAngle(t[e[r]].position));
                        let n = r.meanAngle(i),
                            s = Math.PI / 2;
                        return n = Math.round(Math.round(n / s) * s), 2 === n ? "down" : -2 === n ? "up" : 0 === n || -0 === n ? "right" : 3 === n || -3 === n ? "left" : "down"
                    }
                    getNeighbours(t = null) {
                        if (null === t) return this.neighbours.slice();
                        let e = Array();
                        for (let i = 0; i < this.neighbours.length; i++) this.neighbours[i] !== t && e.push(this.neighbours[i]);
                        return e
                    }
                    getDrawnNeighbours(t) {
                        let e = Array();
                        for (let i = 0; i < this.neighbours.length; i++) t[this.neighbours[i]].value.isDrawn && e.push(this.neighbours[i]);
                        return e
                    }
                    getNeighbourCount() {
                        return this.neighbourCount
                    }
                    getSpanningTreeNeighbours(t = null) {
                        let e = Array();
                        for (let i = 0; i < this.spanningTreeChildren.length; i++) void 0 !== t && t == this.spanningTreeChildren[i] || e.push(this.spanningTreeChildren[i]);
                        return null != this.parentVertexId && (void 0 !== t && t == this.parentVertexId || e.push(this.parentVertexId)), e
                    }
                    getNextInRing(t, e, i) {
                        let r = this.getNeighbours();
                        for (let s = 0; s < r.length; s++)
                            if (n.contains(t[r[s]].value.rings, {
                                value: e
                            }) && r[s] != i) return r[s];
                        return null
                    }
                }
                t.exports = o
            }
        },
        e = {};

        function i(r) {
            if (e[r]) return e[r].exports;
            var n = e[r] = {
                exports: {}
            };
            return t[r](n, n.exports, i), n.exports
        }
        i.n = t => {
            var e = t && t.__esModule ? () => t.default : () => t;
            return i.d(e, {
                a: e
            }), e
        }, i.d = (t, e) => {
            for (var r in e) i.o(e, r) && !i.o(t, r) && Object.defineProperty(t, r, {
                enumerable: !0,
                get: e[r]
            })
        }, i.o = (t, e) => Object.prototype.hasOwnProperty.call(t, e), i.r = t => {
            "undefined" != typeof Symbol && Symbol.toStringTag && Object.defineProperty(t, Symbol.toStringTag, {
                value: "Module"
            }), Object.defineProperty(t, "__esModule", {
                value: !0
            })
        };
        var r = {};
        return (() => {
            "use strict";
            i.r(r), i.d(r, {
                clean2d: () => o
            });
            var t = i(237),
                e = i.n(t),
                n = i(19),
                s = i.n(n);

            function o(t) {
                const i = new (e())({}),
                    r = s().parse(t);
                i.initDraw(r, "light", !1), i.processGraph();
                let n = i.graph.vertices,
                    o = Array();
                for (let t = 0; t < n.length; t++) {
                    let e = n[t].position;
                    o.push([e.x, e.y])
                }
                return o
            }
        })(), r
    })()
}));