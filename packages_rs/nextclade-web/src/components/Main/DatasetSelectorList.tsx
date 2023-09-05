import { isNil } from 'lodash'
import React, { useCallback, useMemo } from 'react'

import { ListGroup, ListGroupItem } from 'reactstrap'
import { useRecoilValue } from 'recoil'
import { minimizerIndexVersionAtom } from 'src/state/dataset.state'
import styled from 'styled-components'

import type { Dataset } from 'src/types'
import { areDatasetsEqual } from 'src/types'
import { search } from 'src/helpers/search'
import { DatasetInfo } from 'src/components/Main/DatasetInfo'

// export const DatasetSelectorContainer = styled.div`
//   flex: 1 0 100%;
//   display: flex;
//   flex-direction: column;
//   overflow: hidden;
//   height: 100%;
//   border: 1px #ccc solid;
//   border-radius: 5px;
// `

const DATASET_AUTODETECT: Dataset = {
  path: 'autodetect',
  enabled: true,
  official: true,
  attributes: {
    name: { value: 'autodetect', valueFriendly: 'Autodetect' },
    reference: { value: 'autodetect', valueFriendly: 'Autodetect' },
  },
  files: {
    reference: '',
    pathogenJson: '',
  },
}

export const DatasetSelectorUl = styled(ListGroup)`
  flex: 1;
  overflow-y: scroll;
  height: 100%;
`

export const DatasetSelectorLi = styled(ListGroupItem)<{ $isDimmed?: boolean }>`
  list-style: none;
  margin: 0;
  padding: 0.5rem;
  cursor: pointer;
  opacity: ${(props) => props.$isDimmed && 0.33};
  background-color: transparent;
`

export interface DatasetSelectorListItemProps {
  dataset: Dataset
  isCurrent?: boolean
  isDimmed?: boolean
  onClick?: () => void
}

export function DatasetSelectorListItem({ dataset, isCurrent, isDimmed, onClick }: DatasetSelectorListItemProps) {
  return (
    <DatasetSelectorLi $isDimmed={isDimmed} aria-current={isCurrent} active={isCurrent} onClick={onClick}>
      <DatasetInfo dataset={dataset} />
    </DatasetSelectorLi>
  )
}

export interface DatasetSelectorListProps {
  datasets: Dataset[]
  searchTerm: string
  datasetHighlighted?: Dataset

  onDatasetHighlighted(dataset?: Dataset): void
}

export function DatasetSelectorList({
  datasets,
  searchTerm,
  datasetHighlighted,
  onDatasetHighlighted,
}: DatasetSelectorListProps) {
  const minimizerIndexVersion = useRecoilValue(minimizerIndexVersionAtom)

  const onItemClick = useCallback((dataset: Dataset) => () => onDatasetHighlighted(dataset), [onDatasetHighlighted])

  const autodetectItem = useMemo(() => {
    if (isNil(minimizerIndexVersion)) {
      return null
    }

    return (
      <DatasetSelectorListItem
        dataset={DATASET_AUTODETECT}
        onClick={onItemClick(DATASET_AUTODETECT)}
        isCurrent={areDatasetsEqual(DATASET_AUTODETECT, datasetHighlighted)}
      />
    )
  }, [datasetHighlighted, minimizerIndexVersion, onItemClick])

  const { itemsStartWith, itemsInclude, itemsNotInclude } = useMemo(() => {
    if (searchTerm.trim().length === 0) {
      return { itemsStartWith: datasets, itemsInclude: [], itemsNotInclude: [] }
    }

    return search(datasets, searchTerm, (dataset) => [
      dataset.attributes.name.value,
      dataset.attributes.name.valueFriendly ?? '',
      dataset.attributes.reference.value,
    ])
  }, [datasets, searchTerm])

  return (
    // <DatasetSelectorContainer>
    <DatasetSelectorUl>
      {autodetectItem}

      {[itemsStartWith, itemsInclude].map((datasets) =>
        datasets.map((dataset) => (
          <DatasetSelectorListItem
            key={dataset.path}
            dataset={dataset}
            onClick={onItemClick(dataset)}
            isCurrent={areDatasetsEqual(dataset, datasetHighlighted)}
          />
        )),
      )}

      {[itemsNotInclude].map((datasets) =>
        datasets.map((dataset) => (
          <DatasetSelectorListItem
            key={dataset.path}
            dataset={dataset}
            onClick={onItemClick(dataset)}
            isCurrent={areDatasetsEqual(dataset, datasetHighlighted)}
            isDimmed
          />
        )),
      )}
    </DatasetSelectorUl>
    // </DatasetSelectorContainer>
  )
}
