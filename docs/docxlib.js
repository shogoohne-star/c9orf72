// Word 文書生成の共通部品
const fs = require('fs');
const path = require('path');
const D = require('docx');
const {
  Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell,
  WidthType, AlignmentType, HeadingLevel, BorderStyle, ShadingType, ExternalHyperlink,
  VerticalAlign, LevelFormat, convertMillimetersToTwip,
} = D;

// ---- 配色（HTML 版のライトテーマから） -------------------------------------
const C = {
  ink: '111C25', ink2: '3F5364', muted: '6E8497',
  rule: 'D9E1E9', rule2: 'C3D0DB', surface2: 'EDF1F5', sunk: 'E7EDF2',
  accent: 'AC2C58', accentSoft: 'FAE8EF',
  steel: '2A6B8C', steelSoft: 'E2EEF5',
  amber: '9E5C28', amberSoft: 'F7EDE3',
};
const JP = { ascii: 'Yu Gothic', eastAsia: 'Yu Gothic', hAnsi: 'Yu Gothic', cs: 'Yu Gothic' };
const MONO = { ascii: 'Consolas', eastAsia: 'MS Gothic', hAnsi: 'Consolas', cs: 'Consolas' };

const CONTENT_W = 9700;          // twip（A4・左右余白 18mm）
const MAX_IMG_PX = 640;

// ---- インライン記法 --------------------------------------------------------
//   **太字**   `等幅`   ^{上付き}   _{下付き}   [文字](URL)
function runs(src, base = {}) {
  if (Array.isArray(src)) return src.flatMap((s) => runs(s, base));
  const out = [];
  const re = /\*\*([^*]+)\*\*|`([^`]+)`|\^\{([^}]*)\}|_\{([^}]*)\}|\[([^\]]+)\]\(([^)]+)\)/g;
  let last = 0, m;
  const push = (text, extra) => {
    if (!text) return;
    out.push(new TextRun({ text, font: JP, ...base, ...extra }));
  };
  while ((m = re.exec(src))) {
    push(src.slice(last, m.index));
    if (m[1] !== undefined) push(m[1], { bold: true });
    else if (m[2] !== undefined) out.push(new TextRun({ text: m[2], font: MONO, size: (base.size || 20) - 2, ...base, font: MONO }));
    else if (m[3] !== undefined) push(m[3], { superScript: true });
    else if (m[4] !== undefined) push(m[4], { subScript: true });
    else out.push(new ExternalHyperlink({
      link: m[6],
      children: [new TextRun({ text: m[5], font: JP, ...base, style: 'Hyperlink' })],
    }));
    last = re.lastIndex;
  }
  push(src.slice(last));
  return out;
}

// ---- ブロック --------------------------------------------------------------
const P = (src, o = {}) => new Paragraph({
  children: runs(src, { size: o.size || 20, color: o.color || C.ink, italics: o.italics }),
  spacing: { before: o.before ?? 0, after: o.after ?? 120, line: o.line ?? 300 },
  alignment: o.align,
  indent: o.indent,
  ...(o.border ? { border: o.border } : {}),
  ...(o.shading ? { shading: o.shading } : {}),
  ...(o.bullet ? { bullet: o.bullet } : {}),
  ...(o.numbering ? { numbering: o.numbering } : {}),
  keepNext: o.keepNext,
});

const H1 = (t) => new Paragraph({
  children: runs(t, { size: 40, bold: true, color: C.ink }),
  spacing: { before: 0, after: 200, line: 340 },
});

const H2 = (t, num) => new Paragraph({
  heading: HeadingLevel.HEADING_1,
  children: [
    ...(num ? runs(num + '   ', { size: 22, bold: true, color: C.accent }) : []),
    ...runs(t, { size: 28, bold: true, color: C.ink }),
  ],
  spacing: { before: 400, after: 160, line: 320 },
  border: { bottom: { style: BorderStyle.SINGLE, size: 12, color: C.ink, space: 6 } },
  keepNext: true,
});

const H3 = (t) => new Paragraph({
  heading: HeadingLevel.HEADING_2,
  children: runs(t, { size: 23, bold: true, color: C.steel }),
  spacing: { before: 280, after: 110, line: 300 },
  keepNext: true,
});

const SUB = (t) => P(t, { size: 18, color: C.muted, after: 160 });

const EYEBROW = (t) => new Paragraph({
  children: runs(t, { size: 17, bold: true, color: C.accent, characterSpacing: 20 }),
  spacing: { after: 100 },
});

const LEDE = (t) => P(t, { size: 21, color: C.ink2, after: 200, line: 330 });

// 引用ボックス（左に色帯 + 淡い背景）
function CALL(title, paras, variant = 'accent') {
  const col = { accent: C.accent, steel: C.steel, amber: C.amber }[variant];
  const bg = { accent: C.accentSoft, steel: C.steelSoft, amber: C.amberSoft }[variant];
  const kids = [];
  if (title) kids.push(P(title, { size: 20, after: 60 }).constructor ? new Paragraph({
    children: runs(title, { size: 20, bold: true, color: col }),
    spacing: { after: 60, line: 300 },
  }) : null);
  paras.forEach((t, i) => kids.push(new Paragraph({
    children: runs(t, { size: 19, color: C.ink }),
    spacing: { after: i === paras.length - 1 ? 0 : 100, line: 300 },
  })));
  return new Table({
    width: { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: [CONTENT_W],
    borders: {
      top: { style: BorderStyle.NONE }, bottom: { style: BorderStyle.NONE },
      right: { style: BorderStyle.NONE }, insideHorizontal: { style: BorderStyle.NONE },
      insideVertical: { style: BorderStyle.NONE },
      left: { style: BorderStyle.SINGLE, size: 18, color: col },
    },
    rows: [new TableRow({
      children: [new TableCell({
        width: { size: CONTENT_W, type: WidthType.DXA },
        shading: { type: ShadingType.CLEAR, fill: bg, color: 'auto' },
        margins: { top: 140, bottom: 140, left: 180, right: 180 },
        children: kids.filter(Boolean),
      })],
    })],
  });
}

const SPACER = (h = 140) => new Paragraph({ text: '', spacing: { after: h, line: 20 } });

// 式ブロック（等幅・淡い背景）
function EQN(lines) {
  const kids = lines.map((l, i) => new Paragraph({
    children: [new TextRun({ text: l, font: MONO, size: 17, color: C.ink })],
    spacing: { after: i === lines.length - 1 ? 0 : 40, line: 260 },
  }));
  return new Table({
    width: { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: [CONTENT_W],
    borders: allBorders(C.rule),
    rows: [new TableRow({
      children: [new TableCell({
        width: { size: CONTENT_W, type: WidthType.DXA },
        shading: { type: ShadingType.CLEAR, fill: C.sunk, color: 'auto' },
        margins: { top: 120, bottom: 120, left: 180, right: 180 },
        children: kids,
      })],
    })],
  });
}

function allBorders(color, size = 4) {
  const b = { style: BorderStyle.SINGLE, size, color };
  return { top: b, bottom: b, left: b, right: b, insideHorizontal: b, insideVertical: b };
}

// ---- 表 --------------------------------------------------------------------
// spec: {caption, head:[...], rows:[[...]], hi:[i], w:[相対幅], n:[数値列index]}
function TABLE(spec) {
  const nCol = spec.head.length;
  const weights = spec.w || new Array(nCol).fill(1);
  const total = weights.reduce((a, b) => a + b, 0);
  let widths = weights.map((x) => Math.round(CONTENT_W * x / total));
  widths[nCol - 1] += CONTENT_W - widths.reduce((a, b) => a + b, 0);
  const nSet = new Set(spec.n || []);
  const fs_ = spec.size || 16;

  const cell = (txt, ci, opts = {}) => new TableCell({
    width: { size: widths[ci], type: WidthType.DXA },
    columnSpan: opts.span,
    shading: opts.fill ? { type: ShadingType.CLEAR, fill: opts.fill, color: 'auto' } : undefined,
    margins: { top: 60, bottom: 60, left: 90, right: 90 },
    verticalAlign: VerticalAlign.TOP,
    children: [new Paragraph({
      children: runs(String(txt), { size: fs_, bold: opts.bold, color: opts.color || C.ink }),
      spacing: { after: 0, line: 250 },
      alignment: nSet.has(ci) && !opts.head ? AlignmentType.RIGHT : AlignmentType.LEFT,
    })],
  });

  const rows = [new TableRow({
    tableHeader: true,
    children: spec.head.map((h, i) => cell(h, i, { head: true, bold: true, fill: C.surface2, color: C.ink2 })),
  })];
  spec.rows.forEach((r, ri) => {
    const hi = (spec.hi || []).includes(ri);
    let ci = 0;
    rows.push(new TableRow({
      children: r.map((v) => {
        const isObj = v && typeof v === 'object';
        const span = isObj ? v.span : 1;
        const c = cell(isObj ? v.t : v, ci, {
          fill: hi ? C.accentSoft : undefined, bold: hi && ci === 0, span,
        });
        ci += span || 1;
        return c;
      }),
    }));
  });

  const out = [];
  if (spec.caption) out.push(new Paragraph({
    children: runs(spec.caption, { size: 17, color: C.muted }),
    spacing: { before: 60, after: 70, line: 260 }, keepNext: true,
  }));
  out.push(new Table({
    width: { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: widths,
    borders: allBorders(C.rule),
    rows,
  }));
  out.push(SPACER(160));
  return out;
}

// ---- 図 --------------------------------------------------------------------
const PNG_DIR = path.join(__dirname, 'png');
const dims = {};
function imgSize(name, maxW = MAX_IMG_PX) {
  if (!dims[name]) {
    const buf = fs.readFileSync(path.join(PNG_DIR, `${name}.png`));
    dims[name] = { w: buf.readUInt32BE(16), h: buf.readUInt32BE(20) };
  }
  const { w, h } = dims[name];
  const s = Math.min(maxW / (w / 3), 1);       // 3 倍解像度で焼いてある
  return { width: Math.round((w / 3) * s), height: Math.round((h / 3) * s) };
}

function FIG(name, caption, legend, maxW) {
  const sz = imgSize(name, maxW);
  const out = [new Paragraph({
    children: [new ImageRun({
      type: 'png',
      data: fs.readFileSync(path.join(PNG_DIR, `${name}.png`)),
      transformation: sz,
    })],
    alignment: AlignmentType.CENTER,
    spacing: { before: 80, after: legend ? 40 : 70 },
    keepNext: true,
  })];
  if (legend) out.push(new Paragraph({
    children: legend.flatMap(([col, lab], i) => [
      new TextRun({ text: (i ? '     ' : '') + '■ ', font: JP, size: 17, color: C[col] }),
      new TextRun({ text: lab, font: JP, size: 17, color: C.ink2 }),
    ]),
    alignment: AlignmentType.CENTER,
    spacing: { after: 70 }, keepNext: true,
  }));
  if (caption) out.push(new Paragraph({
    children: runs(caption, { size: 17, color: C.muted }),
    spacing: { after: 200, line: 270 },
    border: { top: { style: BorderStyle.SINGLE, size: 3, color: C.rule, space: 5 } },
  }));
  out.push(SPACER(60));
  return out;
}

// ---- 箇条書き / 手順 -------------------------------------------------------
const BULLETS = (items) => items.map((t) => new Paragraph({
  children: runs(t, { size: 20, color: C.ink }),
  bullet: { level: 0 },
  spacing: { after: 90, line: 300 },
}));

const STEPS = (items) => items.flatMap(([head, body], i) => [
  new Paragraph({
    children: [
      new TextRun({ text: `${i + 1}. `, font: MONO, size: 20, bold: true, color: C.accent }),
      ...runs(head, { size: 20, bold: true, color: C.ink }),
    ],
    spacing: { before: 120, after: 40, line: 300 }, keepNext: true,
  }),
  new Paragraph({
    children: runs(body, { size: 19, color: C.ink2 }),
    indent: { left: 340 },
    spacing: { after: 90, line: 290 },
  }),
]);

// ---- 用語集（2 列の表） ----------------------------------------------------
function GLOSSARY(items) {
  const W0 = 2400, W1 = CONTENT_W - W0;
  const rows = items.map(([term, def]) => new TableRow({
    children: [
      new TableCell({
        width: { size: W0, type: WidthType.DXA },
        margins: { top: 70, bottom: 70, left: 100, right: 100 },
        children: [new Paragraph({ children: runs(term, { size: 18, bold: true, color: C.ink }), spacing: { after: 0, line: 260 } })],
      }),
      new TableCell({
        width: { size: W1, type: WidthType.DXA },
        margins: { top: 70, bottom: 70, left: 100, right: 100 },
        children: [new Paragraph({ children: runs(def, { size: 18, color: C.ink2 }), spacing: { after: 0, line: 270 } })],
      }),
    ],
  }));
  return [new Table({
    width: { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: [W0, W1],
    borders: allBorders(C.rule),
    rows,
  }), SPACER(160)];
}

// ---- 表紙まわり ------------------------------------------------------------
function METAROW(pairs) {
  const n = pairs.length;
  const w = Math.floor(CONTENT_W / n);
  const widths = new Array(n).fill(w);
  widths[n - 1] += CONTENT_W - w * n;
  return [new Table({
    width: { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: widths,
    borders: {
      top: { style: BorderStyle.SINGLE, size: 4, color: C.rule },
      bottom: { style: BorderStyle.SINGLE, size: 4, color: C.rule },
      left: { style: BorderStyle.NONE }, right: { style: BorderStyle.NONE },
      insideHorizontal: { style: BorderStyle.NONE },
      insideVertical: { style: BorderStyle.SINGLE, size: 4, color: C.rule },
    },
    rows: [new TableRow({
      children: pairs.map(([k, v], i) => new TableCell({
        width: { size: widths[i], type: WidthType.DXA },
        margins: { top: 100, bottom: 100, left: 120, right: 120 },
        children: [
          new Paragraph({ children: runs(k, { size: 15, color: C.muted, bold: true }), spacing: { after: 20, line: 220 } }),
          new Paragraph({ children: runs(v, { size: 17, color: C.ink }), spacing: { after: 0, line: 250 } }),
        ],
      })),
    })],
  }), SPACER(240)];
}

// ---- 要点タイル（2 x N の表） ----------------------------------------------
function STATS(items) {
  const per = 2, W = Math.floor(CONTENT_W / per);
  const rows = [];
  for (let i = 0; i < items.length; i += per) {
    rows.push(new TableRow({
      children: items.slice(i, i + per).map(({ k, v, c, col }) => new TableCell({
        width: { size: W, type: WidthType.DXA },
        margins: { top: 130, bottom: 130, left: 150, right: 150 },
        borders: { top: { style: BorderStyle.SINGLE, size: 18, color: C[col || 'steel'] } },
        children: [
          new Paragraph({ children: runs(k, { size: 16, bold: true, color: C.muted }), spacing: { after: 50, line: 230 } }),
          new Paragraph({ children: runs(v, { size: 32, bold: true, color: C.ink }), spacing: { after: 60, line: 300 } }),
          new Paragraph({ children: runs(c, { size: 17, color: C.ink2 }), spacing: { after: 0, line: 270 } }),
        ],
      })),
    }));
  }
  return [new Table({
    width: { size: CONTENT_W, type: WidthType.DXA },
    columnWidths: new Array(per).fill(W),
    borders: allBorders(C.rule),
    rows,
  }), SPACER(200)];
}

// ---- 構造式カード ----------------------------------------------------------
function STRUCT({ n, title, tag, rows, smiles }) {
  const out = [new Paragraph({
    children: [
      new TextRun({ text: `${n}  `, font: MONO, size: 22, bold: true, color: C.accent }),
      ...runs(title, { size: 23, bold: true, color: C.ink }),
      ...(tag ? runs(`   （${tag}）`, { size: 18, color: C.muted }) : []),
    ],
    spacing: { before: 240, after: 100, line: 300 },
    border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: C.rule2, space: 4 } },
    keepNext: true,
  })];
  rows.forEach(([lab, sub, img]) => {
    out.push(new Paragraph({
      children: [
        ...runs(lab, { size: 19, bold: true, color: C.ink }),
        ...runs('　' + sub, { size: 18, color: C.ink2 }),
      ],
      spacing: { before: 100, after: 60, line: 280 }, keepNext: true,
    }));
    out.push(...FIG(img, null, null, 560));
  });
  smiles.forEach((s) => out.push(new Paragraph({
    children: runs(s, { size: 16, color: C.muted }),
    spacing: { after: 60, line: 260 },
  })));
  out.push(SPACER(120));
  return out;
}

// ---- 文書 ------------------------------------------------------------------
function makeDoc({ title, subject, children }) {
  return new Document({
    creator: 'iedda_template_sim_v18.py',
    title,
    description: subject,
    styles: {
      default: {
        document: { run: { font: JP, size: 20, color: C.ink }, paragraph: { spacing: { line: 300 } } },
        heading1: { run: { font: JP, size: 28, bold: true, color: C.ink } },
        heading2: { run: { font: JP, size: 23, bold: true, color: C.steel } },
      },
      paragraphStyles: [{
        id: 'Normal', name: 'Normal', quickFormat: true,
        run: { font: JP, size: 20, color: C.ink },
        paragraph: { spacing: { line: 300 } },
      }],
    },
    sections: [{
      properties: {
        page: {
          size: { width: convertMillimetersToTwip(210), height: convertMillimetersToTwip(297) },
          margin: {
            top: convertMillimetersToTwip(20), bottom: convertMillimetersToTwip(20),
            left: convertMillimetersToTwip(18), right: convertMillimetersToTwip(18),
            footer: convertMillimetersToTwip(12),
          },
        },
      },
      footers: {
        default: new D.Footer({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [new TextRun({ children: [D.PageNumber.CURRENT], font: JP, size: 16, color: C.muted })],
          })],
        }),
      },
      children,
    }],
  });
}

async function write(doc, file) {
  const buf = await Packer.toBuffer(doc);
  fs.writeFileSync(file, buf);
  console.log(`${file}: ${(buf.length / 1024).toFixed(0)} KB`);
}

module.exports = { C, JP, MONO, CONTENT_W, runs, P, H1, H2, H3, SUB, EYEBROW, LEDE,
  CALL, EQN, TABLE, FIG, BULLETS, STEPS, GLOSSARY, METAROW, STATS, STRUCT, SPACER,
  makeDoc, write, D };
