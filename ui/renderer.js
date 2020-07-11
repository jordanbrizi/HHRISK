const { ipcRenderer } = require('electron')
const path = require('path')
const { count } = require('console')
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
	const arquivos = []
	const jsons = []

	fs.readdirSync(path.resolve('../bin/Results/')).forEach(arquivo => {
		arquivos.push(`${arquivo}`)
	})

	for (j = 0; j < arquivos.length; j++) {
		if (arquivos[j].includes('.json')) jsons.push(arquivos[j])
	}

	return {
		quantidade: () => jsons.length,
		arquivo: chave => jsons[chave],
		pegar: arquivo => require(path.resolve(`../bin/Results/${arquivo}`))
	}
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

module.exports = {pkg, lang, Dados, ipcRenderer}