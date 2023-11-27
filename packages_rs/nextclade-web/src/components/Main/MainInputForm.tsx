import React, { useCallback, useMemo } from 'react'
import { isNil } from 'lodash'
import { useRouter } from 'next/router'
import { useRecoilState, useRecoilValue } from 'recoil'
import { SuggestionAlertMainPage } from 'src/components/Main/SuggestionAlertMainPage'
import styled from 'styled-components'
import { datasetCurrentAtom } from 'src/state/dataset.state'
import { hasRequiredInputsAtom } from 'src/state/inputs.state'
import { shouldSuggestDatasetsOnDatasetPageAtom } from 'src/state/settings.state'
import { useRunSeqAutodetect } from 'src/hooks/useRunSeqAutodetect'
import { useUpdatedDatasetIndex } from 'src/io/fetchDatasets'
import { ButtonChangeDataset, DatasetNoneSection } from 'src/components/Main/ButtonChangeDataset'
import { ButtonRun } from 'src/components/Main/ButtonRun'
import { DatasetCurrentSummary } from 'src/components/Main/DatasetCurrentSummary'
import { MainSectionTitle } from 'src/components/Main/MainSectionTitle'
import { QuerySequenceFilePicker } from 'src/components/Main/QuerySequenceFilePicker'
import { QuerySequenceList } from 'src/components/Main/QuerySequenceList'
import { SelectDatasetHelp } from 'src/components/Help/SelectDatasetHelp'
import { SuggestionPanel } from 'src/components/Main/SuggestionPanel'
import { useRunAnalysis } from 'src/hooks/useRunAnalysis'
import { useTranslationSafe } from 'src/helpers/useTranslationSafe'

const ContainerFixed = styled.div`
  display: flex;
  flex: 1;
  flex-direction: column;
  overflow: hidden;
  width: 100%;
  margin: 0 auto;
  max-width: 1000px;
`

const Container = styled.div`
  display: flex;
  flex: 1;
  flex-direction: column;
  overflow: hidden;
`

const ContainerColumns = styled.div`
  display: flex;
  flex-direction: row;
  overflow: hidden;
`

const Header = styled.div`
  display: flex;
  flex: 0;
  padding-left: 10px;
  margin-top: 10px;
  margin-bottom: 3px;
`

const Main = styled.div`
  display: flex;
  flex-direction: column;
  overflow: hidden;
`

const Footer = styled.div`
  display: flex;
  flex: 1;
`

const FooterWithoutOverflow = styled(Footer)`
  overflow: hidden;
`

export function Landing() {
  // This periodically fetches dataset index and updates the list of datasets.
  useUpdatedDatasetIndex()

  const { push } = useRouter()
  const runAutodetect = useRunSeqAutodetect()
  const hasRequiredInputs = useRecoilValue(hasRequiredInputsAtom)
  const shouldSuggestDatasets = useRecoilValue(shouldSuggestDatasetsOnDatasetPageAtom)

  const toDatasetSelection = useCallback(() => {
    void push('/dataset') // eslint-disable-line no-void
    if (shouldSuggestDatasets && hasRequiredInputs) {
      runAutodetect()
    }
  }, [hasRequiredInputs, push, runAutodetect, shouldSuggestDatasets])

  return (
    <ContainerFixed>
      <Header>
        <MainSectionTitle />
      </Header>

      <Main className="mt-4">
        <ContainerColumns>
          <QuerySequenceFilePicker />
          <DatasetCurrentOrSelectButton toDatasetSelection={toDatasetSelection} />
        </ContainerColumns>
      </Main>

      <FooterWithoutOverflow>
        <Container>
          <Main>
            <QuerySequenceList />
          </Main>
        </Container>
      </FooterWithoutOverflow>
    </ContainerFixed>
  )
}

export interface DatasetCurrentOrSelectButtonProps {
  toDatasetSelection(): void
}

function DatasetCurrentOrSelectButton({ toDatasetSelection }: DatasetCurrentOrSelectButtonProps) {
  const { t } = useTranslationSafe()
  const run = useRunAnalysis()

  const [dataset, _0] = useRecoilState(datasetCurrentAtom)

  const text = useMemo(() => {
    if (isNil(dataset)) {
      return t('Select dataset')
    }
    return t('Selected dataset')
  }, [dataset, t])

  if (!dataset) {
    return (
      <Container>
        <Header>
          <Title>
            <H4Inline>{text}</H4Inline>
            <SelectDatasetHelp />
          </Title>
        </Header>

        <Main>
          <DatasetNoneSection toDatasetSelection={toDatasetSelection} />
        </Main>

        <Footer>
          <div>
            <SuggestionPanel />
            <SuggestionAlertMainPage />
          </div>
        </Footer>
      </Container>
    )
  }

  return (
    <Container>
      <Header>
        <Title>
          <H4Inline>{text}</H4Inline>
          <SelectDatasetHelp />
        </Title>
      </Header>

      <Main>
        <DatasetCurrentSummary />
      </Main>

      <Footer>
        <div className="w-100 d-flex flex-column">
          <SuggestionPanel />
          <SuggestionAlertMainPage />
        </div>
      </Footer>

      <Footer>
        <ButtonChangeDataset className="mr-auto my-2" onClick={toDatasetSelection} />
        <ButtonRun className="ml-auto my-2" onClick={run} />
      </Footer>
    </Container>
  )
}

const Title = styled.h4`
  display: flex;
  flex: 1;
`

const H4Inline = styled.h4`
  display: inline-flex;
  margin: auto 0;
`
