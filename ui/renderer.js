const { ipcRenderer } = require('electron')
const { createBrotliDecompress } = require('zlib')
const { cachedDataVersionTag } = require('v8')
const defaultLang = Intl.DateTimeFormat().resolvedOptions().locale
const pkg = () => require('../package')

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
	const path = require('path')
	const fs = require('fs')
	const files = []
	fs.readdirSync(path.resolve('./bin/Results/')).forEach(arquivo => {
		files.push(`${arquivo}`)
	})

	const jsons = files.filter(a => a.includes('.json'))

	const hhr = ['Doses, HQ and CR.json',
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
		pegar: arquivo => require(path.resolve(`./bin/Results/${arquivo}`))
	}
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

const gerarOds = t => {
	const path = require('path')
	const xlsx = require('xlsx')
	const caminho = path.resolve(`./bin/Results/`)
	const planilhas = []
	jsons = Dados().arquivos()[t] // PEGA OS JSONS DO GRUPO t
	jsons.forEach(json => {
		arquivo = Dados().pegar(json)
		chaves = Object.keys(arquivo)
		const wb = xlsx.utils.book_new()
		chaves.forEach(chave => {
			chaveNew = chave.substring(0, 28) + '...'
			keys = Object.keys(arquivo[chave][0])
			ws = xlsx.utils.json_to_sheet(arquivo[chave])
			//const merge = [{ s: { r: 0, c: 0 }, e: { r: 0, c: (keys.length -1) } }]
			//ws["!merges"] = merge
			xlsx.utils.book_append_sheet(
				wb,
				ws,
				chaveNew
			)
		})
		let filePath = `./bin/Results/${json.replace('.json', '')}.ods`
		xlsx.writeFile(wb, filePath)
		planilhas.push(path.resolve(filePath))
	})
	
	ipcRenderer.send('salvar_planilha', {path: caminho, sheets: planilhas})
}
module.exports = {pkg, lang, Dados, ipcRenderer, gerarOds}