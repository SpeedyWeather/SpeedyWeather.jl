import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import path from 'path'

const mathjax = mathjaxPlugin()

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const navTemp: { nav: any } = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker',
  }
]

// Build a per-section sidebar: each page maps to only its own section's items
function buildSidebar(nav: any[]): Record<string, any[]> {
  const result: Record<string, any[]> = {}
  for (const section of nav) {
    if (!section.items) continue  // skip flat items like Home, API, VersionPicker
    const sectionSidebar = [{ text: section.text, collapsed: false, items: section.items }]
    for (const page of section.items) {
      if (page.link) {
        result[page.link] = sectionSidebar
      }
    }
  }
  return result
}

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS',

  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', { src: '/SpeedyWeatherDocumentation/versions.js' }],
    ['script', { src: `${baseTemp.base}siteinfo.js` }],
  ],

  ignoreDeadLinks: true,

  markdown: {
    config(md) {
      md.use(tabsMarkdownPlugin);
      mathjax.markdownConfig(md);
    },
    theme: {
      light: 'github-light',
      dark: 'github-dark',
    },
  },

  vite: {
    plugins: [
      mathjax.vitePlugin,
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components'),
      },
    },
    build: {
      assetsInlineLimit: 0,
    },
    optimizeDeps: {
      exclude: [
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ],
    },
    ssr: {
      noExternal: [
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ],
    },
  },

  themeConfig: {
    outline: 'deep',
    logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    search: {
      provider: 'local',
      options: {
        detailedView: true,
      },
    },
    nav,
    sidebar: buildSidebar(navTemp.nav),
    editLink: {
      pattern: 'https://github.com/SpeedyWeather/SpeedyWeather.jl/edit/main/docs/src/:path',
    },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/SpeedyWeather/SpeedyWeather.jl' },
    ],
    footer: {
      message: 'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> and <a href="https://luxdl.github.io/DocumenterVitepress.jl/stable" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>Released under the MIT License. Powered by the <a href="https://www.julialang.org">Julia Programming Language</a>.<br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`,
    },
  },
})
