export function setupKatex() {
    renderMathInElement(document.body, {
        delimiters: [
        {left: "\\begin{equation}", right: "\\end{equation}", display: true},
        {left: "\\begin{align}", right: "\\end{align}", display: true},
        {left: "$$", right: "$$", display: true},
        {left: "$", right: "$", display: false},
        {left: "\\[", right: "\\]", display:true},
        {left: "\\(", right: "\\)", display:false}
        ],
        macros: {
            "\\F": "\\mathbb{F}",
            "\\ps": "\\left\\langle #1 , #2 \\right\\rangle",
            "\\K": "\\mathbb{K}",
            "\\LL": "\\mathscr{L}",
            "\\nint": "\\left\\lceil #1 \\right\\rfloor",
            "\\ps": "\\left\\langle #1 , #2 \\right\\rangle",
            "\\norm": "\\left\\| #1 \\right\\|",
            "\\gscoeff": "\\frac{\\ps{#1}{#2}}{\\norm{#2}^2} #2",
            "\\OO": "\\mathcal{O}",
            "\\ra": "\\rightarrow"
        }
    });
}


/* il serait temps que j y bosse la dessus*/