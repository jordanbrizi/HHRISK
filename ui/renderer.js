const { ipcRenderer } = require('electron')
const path = require('path')
const defaultLang = Intl.DateTimeFormat().resolvedOptions().locale
const pkg = () => require('./package')

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

const lang = (attr = defaultLang) => {
	if(attr != null) {
		return require(`./lang/${attr}`)
	}
	else {
		return require(`./lang/${defaultLang}`)	
	}
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

const Dados = () => {
	const fs = require('fs')
	const files = []

	fs.readdirSync(path.resolve('../bin/Results/')).forEach(arquivo => {
		files.push(`${arquivo}`)
	})

	const jsons = files.filter(a => a.includes('.json'))

	const hhr = ['Doses, HQ and CR.json',
		'Summary Aggregated Risk.json',
		'Extensive Aggregated Risk.json',
		'Summary Cumulative Risk.json',
		'Extensive Cumulative Risk.json',
		'Complementary Analyzes.json'
	] // Human Health Risk
	const er = ['Combined.json', 'Individual.json'] // Ecological Risk
	const rr = ['Radiological Risk.json'] // Radiological Risk

	const arquivos = {
		arquivos_hhr: jsons.filter(a => hhr.includes(a)),
		arquivos_er: jsons.filter(a => er.includes(a)),
		arquivos_rr: jsons.filter(a => rr.includes(a))
	}

	return {
		quantidade: () => jsons.length,
		arquivos: () => jsons,
		pegar: arquivo => require(path.resolve(`../bin/Results/${arquivo}`))
	}
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

module.exports = {pkg, lang, Dados, ipcRenderer}