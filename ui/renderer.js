const { ipcRenderer, remote, shell } = require('electron')
const { createBrotliDecompress } = require('zlib')
const { cachedDataVersionTag } = require('v8')
const defaultLang = Intl.DateTimeFormat().resolvedOptions().locale
const pkg = () => require('../package')
const resultsPath = arg => remote.app.getAppPath() + '\\bin\\Results\\' + arg
const mainPath = arg => remote.app.getAppPath() + '\\' + arg

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

const lang = (attr = defaultLang) => {
	if(attr != null) {
		return require(mainPath(`ui\\lang\\${attr}`))
	}
	else {
		return require(mainPath(`ui\\lang\\${defaultLang}`))	
	}
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

const Dados = () => {
	const path = require('path')
	const fs = require('fs')
	const files = []
	fs.readdirSync(mainPath('\\bin\\Results')).forEach(arquivo => {
		files.push(`${arquivo}`)
	})

	const jsons = files.filter(a => a.includes('.json'))

	const hhr = [
		'Summary Aggregate Risk.json',
		'Extensive Aggregate Risk.json',
		'Summary Cumulative Risk.json',
		'Extensive Cumulative Risk.json',
		'Complementary Analyzes.json'
	] // Human Health Risk
	const er = ['Combined.json', 'Individual.json'] // Ecological Risk
	const rr = ['Radiological Risk.json'] // Radiological Risk

	const arquivos = {
		hhr: jsons.filter(a => hhr.includes(a)),
		er: jsons.filter(a => er.includes(a)),
		rr: jsons.filter(a => rr.includes(a))
	}

	return {
		quantidade: () => arquivos.length,
		arquivos: () => arquivos,
		pegar: arquivo => require(resultsPath(arquivo))
	}
}

module.exports = {pkg, lang, Dados, ipcRenderer, shell, remote}